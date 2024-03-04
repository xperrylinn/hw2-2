#include "common.h"
#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <unordered_set>

/*
* Static global variables
*/
struct grid_cell_t {
    int rank;
    int index;
    std::vector<particle_t*> particles;
    std::vector<particle_t*> ghost_particles;
    std::vector<int> is_ghost_cell_to; // stores rank indices
    
    grid_cell_t* neighbor_up;
    grid_cell_t* neighbor_down;    
    grid_cell_t* neighbor_left_up;
    grid_cell_t* neighbor_right_up;
    grid_cell_t* neighbor_left;
    grid_cell_t* neighbor_right;
    grid_cell_t* neighbor_left_down;
    grid_cell_t* neighbor_right_down;

    // Constructor to initialize member vectors
    grid_cell_t() 
        : rank(0), 
          index(0),
          particles(),  // Initialize particle vector (empty)
          ghost_particles(),  // Initialize ghost_particle vector (empty)
          is_ghost_cell_to(), // Initialize is_ghost_cell_to vector (empty)
          neighbor_up(nullptr),
          neighbor_down(nullptr),
          neighbor_left_up(nullptr),
          neighbor_right_up(nullptr),
          neighbor_left(nullptr),
          neighbor_right(nullptr),
          neighbor_left_down(nullptr),
          neighbor_right_down(nullptr) {}

};

static int grid_dimension;
static double grid_cell_length;
static std::vector<grid_cell_t*> grid_cells;
static std::vector<grid_cell_t*> rank_grid_cells;

/*
* Helper Functions
*/

// Helper function for returning the index of the block that a particle belongs to.
int get_block_index(particle_t* p) {
    return floor(p->x / cutoff) + floor(p->y / cutoff) * grid_dimension;
}

// Helper function for returning the rank of that a grid cell belongs.
int get_grid_cell_rank(int grid_index) {
    return 0;
}

// Helper function for returning the rank of that a grid cell belongs.
int get_rank_from_particle_position(int grid_index) {
    return 0;
}

/*
*   Core simulation functions
*/

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;

    //TODO: 2-way update
    // neighbor.ax -= coef * dx;
    // neighbor.ay -= coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}

void init_simulation(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
	// You can use this space to initialize data objects that you may need
	// This function will be called once before the algorithm begins
	// Do not do any particle simulation here

    // Setup grids each with cutoff-by-cutoff dimensions.
    grid_dimension = ceil(size / cutoff);
    int total_grid_count = grid_dimension * grid_dimension;
    grid_cell_length = size / grid_dimension;

    // Create grid of grid_cell_t
    grid_cells.resize(total_grid_count);
    for (int i = 0; i < total_grid_count; ++i) {
        // Create a new grid_cell_t instance and assign its pointer to the vector
        grid_cells[i] = new grid_cell_t();
    }
    

    // Iterate over all pariticles and place them into their corresponding grid_cell_t
    for (int i = 0; i < num_parts; i++) {
        particle_t* particle = parts + i;
        int block_index = get_block_index(particle);
        fflush(stdout);
        grid_cells[block_index]->particles.push_back(particle);
    }

    // Link neighboring grid_cell_t
    for (int grid_index = 0; grid_index < total_grid_count; grid_index++) {
        int col_index = grid_index % grid_dimension;
        int row_index = grid_index / grid_dimension;
        // Up
        if (row_index > 0) {
            int up_grid_index = grid_index - grid_dimension;
            grid_cells[grid_index]->neighbor_up = grid_cells[up_grid_index];
        }
        // Down
        if (row_index + 1 < grid_dimension) {
            int down_grid_index = grid_index + grid_dimension;
            grid_cells[grid_index]->neighbor_down = grid_cells[down_grid_index];
        }
        // Left
        if (col_index - 1 >= 0) {
            int left_grid_index = grid_index - 1;
            grid_cells[grid_index]->neighbor_left = grid_cells[left_grid_index];
        }
        // Right
        if (col_index + 1 < grid_dimension) {
            int right_grid_index = grid_index + 1;
            grid_cells[grid_index]->neighbor_right = grid_cells[right_grid_index];
        }
        // Up-Left
        if (col_index - 1 >= 0 && row_index > 0) {
            int up_left_grid_index = grid_index - grid_dimension - 1;
            grid_cells[grid_index]->neighbor_left_up = grid_cells[up_left_grid_index];
        }
        // Up-Right
        if (col_index + 1 < grid_dimension && row_index > 0) {
            int up_right_grid_index = grid_index - grid_dimension + 1;
            grid_cells[grid_index]->neighbor_right_up = grid_cells[up_right_grid_index];
        }
        // Down-Left
        if (col_index - 1 >= 0 && row_index + 1 < grid_dimension) {
            int down_left_grid_index = grid_index + grid_dimension - 1;
            grid_cells[grid_index]->neighbor_left_down = grid_cells[down_left_grid_index];
        }
        // Down-Right
        if (col_index + 1 < grid_dimension && row_index + 1 < grid_dimension) {
            int up_right_grid_index = grid_index + grid_dimension + 1;
            grid_cells[grid_index]->neighbor_right_down = grid_cells[up_right_grid_index];
        }
    }

    // Assign grid_cell_t for each rank
    for (int k = 0; k < num_procs; k += 1) {
        for (int i = k * grid_dimension; i < grid_dimension * grid_dimension; i += grid_dimension * num_procs) {
            for (int j = 0; (j < grid_dimension) && (i + j < grid_dimension * grid_dimension); j += 1) {
                // Logging for debugging
                // printf("index: %d, rank: %d, grid_dimension: %d, num_procs: %d, num_parts: %d,\n", i + j, rank, grid_dimension, num_procs, num_parts);
                // fflush(stdout);

                // Compute current grid index and get the corresponding grid_cell_t
                int curr_grid_index = i + j;
                grid_cell_t* curr_grid_cell = grid_cells[curr_grid_index];

                // Assign rank and index to current grid_cell
                curr_grid_cell->rank = k;
                curr_grid_cell->index = curr_grid_index;

                // Add grid_cell to this rank's vector of grid_cell_t
                if (rank == k) {
                    rank_grid_cells.push_back(curr_grid_cell);
                }
            }
        }
    }


    // TODO: Assign ghost particles
    for (int i = 0; i < rank_grid_cells.size(); i += 1) {
        grid_cell_t* curr_grid_cell = rank_grid_cells[i];
        // Up
        if (curr_grid_cell->neighbor_up != NULL && curr_grid_cell->neighbor_up->rank != rank) {
            curr_grid_cell->ghost_particles.insert(
                curr_grid_cell->ghost_particles.end(),
                curr_grid_cell->neighbor_up->particles.begin(),
                curr_grid_cell->neighbor_up->particles.end()
            );
            curr_grid_cell->is_ghost_cell_to.push_back(curr_grid_cell->neighbor_up->index);
        }
        // Down
        if (curr_grid_cell->neighbor_down != NULL && curr_grid_cell->neighbor_down->rank != rank) {
            curr_grid_cell->ghost_particles.insert(
                curr_grid_cell->ghost_particles.end(), 
                curr_grid_cell->neighbor_down->particles.begin(), 
                curr_grid_cell->neighbor_down->particles.end()
            );
            curr_grid_cell->is_ghost_cell_to.push_back(curr_grid_cell->neighbor_down->index);
        }
        // Left
        if (curr_grid_cell->neighbor_left != NULL && curr_grid_cell->neighbor_left->rank != rank) {
            curr_grid_cell->ghost_particles.insert(
                curr_grid_cell->ghost_particles.end(),
                curr_grid_cell->neighbor_left->particles.begin(),
                curr_grid_cell->neighbor_left->particles.end()
            );
            curr_grid_cell->is_ghost_cell_to.push_back(curr_grid_cell->neighbor_left->index);
        }
        // Right
        if (curr_grid_cell->neighbor_right != NULL && curr_grid_cell->neighbor_right->rank != rank) {
            curr_grid_cell->ghost_particles.insert(
                curr_grid_cell->ghost_particles.end(),
                curr_grid_cell->neighbor_right->particles.begin(),
                curr_grid_cell->neighbor_right->particles.end()
            );
            curr_grid_cell->is_ghost_cell_to.push_back(curr_grid_cell->neighbor_right->index);  
        }
        // Up-Left
        if (curr_grid_cell->neighbor_left_up != NULL && curr_grid_cell->neighbor_left_up->rank != rank) {
            curr_grid_cell->ghost_particles.insert(
                curr_grid_cell->ghost_particles.end(),
                curr_grid_cell->neighbor_left_up->particles.begin(),
                curr_grid_cell->neighbor_left_up->particles.end()
            );
            curr_grid_cell->is_ghost_cell_to.push_back(curr_grid_cell->neighbor_left_up->index);
        }
        // Up-Right
        if (curr_grid_cell->neighbor_right_up != NULL && curr_grid_cell->neighbor_right_up->rank != rank) {
            curr_grid_cell->ghost_particles.insert(
                curr_grid_cell->ghost_particles.end(),
                curr_grid_cell->neighbor_right_up->particles.begin(),
                curr_grid_cell->neighbor_right_up->particles.end()
            );
            curr_grid_cell->is_ghost_cell_to.push_back(curr_grid_cell->neighbor_right_up->index);
        }
        // Down-Left
        if (curr_grid_cell->neighbor_left_down != NULL && curr_grid_cell->neighbor_left_down->rank != rank) {
            curr_grid_cell->ghost_particles.insert(
                curr_grid_cell->ghost_particles.end(),
                curr_grid_cell->neighbor_left_down->particles.begin(),
                curr_grid_cell->neighbor_left_down->particles.end()
            );
            curr_grid_cell->is_ghost_cell_to.push_back(curr_grid_cell->neighbor_left_down->index);
        }
        // Down-Right
        if (curr_grid_cell->neighbor_right_down != NULL && curr_grid_cell->neighbor_right_down->rank != rank) {
            curr_grid_cell->ghost_particles.insert(
                curr_grid_cell->ghost_particles.end(), 
                curr_grid_cell->neighbor_right_down->particles.begin(),
                curr_grid_cell->neighbor_right_down->particles.end()
            );
            curr_grid_cell->is_ghost_cell_to.push_back(curr_grid_cell->neighbor_right_down->index);
        }
    }
    
    // Log the part ids owned by each rank for debugging
    // printf("rank %d has the following parts: ", rank);
    // fflush(stdout);
    // for (int element : rank_part_ids) {
    //     printf("%d ", element);
    //     fflush(stdout);
    // }
    // printf("\n");
    // fflush(stdout);

    // More logging to confirm end of fxn was reached
    printf("made it to end of init_simulation\n");
}

void simulate_one_step(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    // Write this function
    printf("simulate_one_step, rank: %d\n", rank);
    fflush(stdout);

    // For each grid cell object calculate forces between all particles
        // For particles within current grid cell
            // Intercel force 
                // Apply force
            // Intracel force
                // For neighbor particle
                    // Apply force

    // For each grid cell object
        // Check if ghost region of other rank
            // If True across rank communication
        // For each particle
            // Move
              // Update particle P to respective grid 
              	// Some particles will no longer be here, we can update regardless or do a check, doesn't matter
            // Afrer move() there will be communication with diff rank
              // In implementation, tally up changes (amount + particle_ids for each neighbor rank)
            // 1. (Universal) If particle moving into another chunk (owned by another rank)
            // 2. (check start grid) If particle started in a gridcell that is in ghost regions of other rank(s)
            // 3. (check end grid) If particle moved into a gridcell that is in ghost regions of other rank(s)
            // Note: two checks for ghost regions, one for starting grid and one for ending grid.
  	
  	// Communication section
  		// For each neighbor rank
  			// Put particles meant for that rank into buffer, send to that rank
  			// MPI_Probe the amount of particles to receive, allocate buffer, receive particles
  			// update local particles in PARTS
  
  	// Update rank-local grid_cells
  		// Clear the previous 
  		// Loop through particles received

}

void gather_for_save(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    // Write this function such that at the end of it, the master (rank == 0)
    // processor has an in-order view of all particles. That is, the array
    // parts is complete and sorted by particle id.

    // Log how many particles were passerd in for debugging
    // for (int i = 0; i < num_parts; i += 1) {
    //     printf("id: %d\n", parts[i].id);
    // }
    printf("made it to beggining of gather_for_save\n");
    fflush(stdout);

    // Vectors for gathering particles from processors
    std::vector<particle_t*> sending_parts;  // Size depends number of particles owned by each processor
    std::vector<particle_t*> receiving_parts(num_parts); // Size

    // Add particles from processor to be sent for gathering
    for (int i = 0; i < rank_grid_cells.size(); i += 1) {
        std::vector<particle_t*> curr_grid_cell_parts = rank_grid_cells[i]->particles;
        for (int j = 0; j < curr_grid_cell_parts.size(); j += 1) {
            sending_parts.push_back(curr_grid_cell_parts[j]);
        }
    }

    // Use variable gather due to particle count varying for each processor
    // Create arrays for specifying count and offsets of particles for each processor
    int* receiving_counts;
    int* offsets;
    if (rank == 0) {
        receiving_counts = new int[num_procs];
        offsets = new int[num_procs];
    }

    // Use gather to populate array for specifying count of partilces for each offset
    printf("rank %d before MPI_Gather\n", rank);
    fflush(stdout);
    int sending_count = sending_parts.size();
    MPI_Gather(&sending_count, 1, MPI_INT, receiving_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Build up offsets array
    if (rank == 0) {
        offsets[0] = 0;
        for (int i = 1; i < num_parts; i += 1) {
            offsets[i] = offsets[i - 1] + receiving_counts[i - 1];
        }
    }

    // Variable gather of particles from processors
    printf("rank %d before MPI_Gatherv\n", rank);
    fflush(stdout);
    MPI_Gatherv(sending_parts.data(), sending_parts.size(), PARTICLE, receiving_parts.data(), receiving_counts, offsets, PARTICLE, 0, MPI_COMM_WORLD);

    // Create in-order view of all particles, sorted by particle id
    if (rank == 0) {
        for (int i = 0; i < num_parts; i += 1) {
            particle_t* curr_part = receiving_parts[i];
            parts[curr_part->id - 1] = *curr_part;
        }
    }

    printf("made it to end of gather_for_save\n");
    fflush(stdout);
}