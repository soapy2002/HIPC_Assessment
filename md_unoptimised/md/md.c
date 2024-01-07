#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "args.h"
#include "boundary.h"
#include "data.h"
#include "setup.h"
#include "vtk.h"

/**
 * @brief This routine calculates the acceleration felt by each particle based on evaluating the Lennard-Jones 
 *        potential with its neighbours. It only evaluates particles within a cut-off radius, and uses cells to 
 *        reduce the search space. It also calculates the potential energy of the system. 
 * 
 * @return double The potential energy
 */
double comp_accel() {
	// zero acceleration for every particle
	for (int i = 1; i < x+1; i++) {
		for (int j = 1; j < y+1; j++) {
			struct particle_t * p = cells[i][j].head;
			while (p != NULL) {
				p->ax = 0.0;
				p->ay = 0.0;
				p = p->next;
			}
		}
	}

	double pot_energy = 0.0;

	for (int i = 1; i < x+1; i++) {
		for (int j = 1; j < y+1; j++) {
			double cell_offset_x = (i-1) * cell_size;
			double cell_offset_y = (j-1) * cell_size;
			struct particle_t * p = cells[i][j].head;
			while (p != NULL) {
				// Compare each particle with all particles in the 9 cells
				for (int a = -1; a <= 1; a++) {
					for (int b = -1; b <= 1; b++) {
						struct particle_t * q = cells[i+a][j+b].head;
						while (q != NULL) {
							// if p and q are the same particle, skip
							if (p == q) {
								q = q->next;
								continue;
							}

							// since particles are stored relative to their cell, calculate the
							// actual x and y coordinates.
							double p_real_x = (cell_offset_x) + p->x;
							double p_real_y = (cell_offset_y) + p->y;
							double q_real_x = ((i+a-1) * cell_size) + q->x;
							double q_real_y = ((j+b-1) * cell_size) + q->y;
							
							// calculate distance in x and y, then absolute distance
							double dx = p_real_x - q_real_x;
							double dy = p_real_y - q_real_y;
							double r_2 = dx*dx + dy*dy;
							
							// if distance less than cut off, calculate force and 
							// use this to calculate acceleration in each dimension
							// calculate potential energy of each particle at the same time
							if (r_2 < r_cut_off_2) {
								double r_2_inv = 1.0 / r_2;
								double r_6_inv = r_2_inv * r_2_inv * r_2_inv;
								
								double f = (48.0 * r_2_inv * r_6_inv * (r_6_inv - 0.5));

								p->ax += f*dx;
								p->ay += f*dy;

								pot_energy += 4.0 * r_6_inv * (r_6_inv - 1.0) - Uc - Duc * (sqrt(r_2) - r_cut_off);
							}
							q = q->next;
						}
					}
				}
				p = p->next;
			}
		}
	}
	// return the average potential energy (i.e. sum / number)
	return pot_energy / num_particles;
}

/**
 * @brief This routine updates the velocity of each particle for half a time step and then 
 *        moves the particle for a whole time step
 * 
 */
void move_particles() {
	// move all particles half a time step
	for (int i = 1; i < x+1; i++) {
		for (int j = 1; j < y+1; j++) {
			struct particle_t * p = cells[i][j].head;
			while (p != NULL) {
				// update velocity to obtain v(t + Dt/2)
				p->vx += dth * p->ax;
				p->vy += dth * p->ay;

				// update particle coordinates to p(t + Dt) (scaled to the cell_size)
				p->x += (dt * p->vx);
				p->y += (dt * p->vy);

				p = p->next;
			}
		}
	}
}

/**
 * @brief This routine updates the cell lists. If a particles coordinates are not within a cell
 *        any more, this function calculates the cell it should be in and performs the move.
 *        If a particle moves more than 1 cell in any direction, this indicates poor settings
 *        and therefore an error is generated.
 * 
 */
void update_cells() {
	// move particles that need to move cell lists
	for (int i = 1; i < x+1; i++) {
		for (int j = 1; j < y+1; j++) {
			struct particle_t * p = cells[i][j].head;
			struct particle_t * p_next;
			while (p != NULL) {
				// we have to store the next particle here, as the remove/add at the end may be destructive
				p_next = p->next;
				// if a particles x or y value is greater than the cell size or less than 0, it must have moved cell
				// do a quick check to make sure its not moved 2 cells (since this means our time step is too large, or something else is going wrong)
				if ((p->x < 0.0) | (p->x >= cell_size) | (p->y < 0.0) | (p->y >= cell_size)) {
					if ((p->x < (-cell_size)) || (p->x >= (2*cell_size)) || (p->y < (-cell_size)) || (p->y >= (2*cell_size))) {
						fprintf(stderr, "A particle has moved more than one cell!\n");
						exit(1);
					}

					// work out whether we've moved a cell in the x and the y dimension
					int x_shift = (p->x < 0.0) ? -1 : (p->x >= cell_size) ? +1 : 0;
					int y_shift = (p->y < 0.0) ? -1 : (p->y >= cell_size) ? +1 : 0;
					
					// the new i and j are +/- 1 in each dimension,
					// but if that means we go out of simulation bounds, wrap it to x and 1
					int new_i = i+x_shift;
					if (new_i == 0) { new_i = x; }
					if (new_i == x+1) { new_i = 1; }
					int new_j = j+y_shift;
					if (new_j == 0) { new_j = y; }
					if (new_j == y+1) { new_j = 1; }
					// update x and y coordinates (i.e. remove the additional cell size)
					p->x = p->x + (x_shift * -cell_size);
					p->y = p->y + (y_shift * -cell_size);

					// remove the particle from its current cell list, then add it to the new cell list
					remove_particle(&(cells[i][j]), p);
					add_particle(&(cells[new_i][new_j]), p);
				}
				p = p_next;
			}
		}
	}
}

/**
 * @brief This updates the velocity of particles for the whole time step (i.e. adds the acceleration for another
 *        half step, since its already done half a time step in the move_particles routine). Additionally, this
 *        function calculated the kinetic energy of the system.
 * 
 * @return double The kinetic energy
 */
double update_velocity() {
	double kinetic_energy = 0.0;

	for (int i = 1; i < x+1; i++) {
		for (int j = 1; j < y+1; j++) {
			struct particle_t * p = cells[i][j].head;
			while (p != NULL) {
				// update velocity again by half time to obtain v(t + Dt)
				p->vx += dth * p->ax;
				p->vy += dth * p->ay;

				// calculate the kinetic energy by adding up the squares of the velocities in each dim
				kinetic_energy += (p->vx * p->vx) + (p->vy * p->vy);

				p = p->next;
			}
		}
	}

	// KE = (1/2)mv^2
	kinetic_energy *= (0.5 / num_particles);
	return kinetic_energy;
}

/**
 * @brief This is the main routine that sets up the problem space and then drives the solving routines.
 * 
 * @param argc The number of arguments passed to the program
 * @param argv An array of the arguments passed to the program
 * @return int The exit code of the application
 */
int main(int argc, char *argv[]) {
	clock_t start_time, end_time;

	start_time = clock();

	// Set default parameters
	set_defaults();
	// parse the arguments
	parse_args(argc, argv);
	// call set up to update defaults
	setup();

	if (verbose) print_opts();
	
	// set up problem
	problem_setup();

	// apply boundary condition (i.e. update pointers on the boundarys to loop periodically)
	apply_boundary();
	
	comp_accel();

	double potential_energy = 0.0;
	double kinetic_energy = 0.0;

	double placeholder = 2.0 / 3.0;

	int iters = 0;
	double t;
	for (t = 0.0; t < t_end; t+=dt, iters++) {
		// move particles half a time step
		move_particles();

		// update cell lists (i.e. move any particles between cell lists if required)
		update_cells();

		// update pointers (because the previous operation might break boundary cell lists)
		apply_boundary();
		
		// compute acceleration for each particle and calculate potential energy
		potential_energy = comp_accel();

		// update velocity based on the acceleration and calculate the kinetic energy
		kinetic_energy = update_velocity();
	
		if (iters % output_freq == 0) {
			// calculate temperature and total energy
			double total_energy = kinetic_energy + potential_energy;
			double temp = kinetic_energy * placeholder;

			printf("Step %8d, Time: %14.8e (dt: %14.8e), Total energy: %14.8e (p:%14.8e,k:%14.8e), Temp: %14.8e\n", iters, t+dt, dt, total_energy, potential_energy, kinetic_energy, temp);
 
			// if output is enabled and checkpointing is enabled, write out
            if ((!no_output) && (enable_checkpoints))
                write_checkpoint(iters, t+dt);
		}
	}

	// calculate the final energy and write out a final status message
	double final_energy = kinetic_energy + potential_energy;
	printf("Step %8d, Time: %14.8e, Final energy: %14.8e\n", iters, t, final_energy);
    printf("Simulation complete.\n");

	// if output is enabled, write the mesh file and the final state
	if (!no_output) {
		write_mesh();
		write_result(iters, t);
	}

	end_time = clock();

	double total_cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

	printf("total cpu time: %lf seconds \n", total_cpu_time);

	return 0;
}

