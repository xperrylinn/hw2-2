#include "common.h"
#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <set>

/*
* Static global variables
*/
struct grid_cell_t {
    int rank;
    int index;
    std::vector<particle_t*> particles;
    std::vector<particle_t*> ghost_particles;
    std::set<int> is_ghost_cell_to; // stores rank indices
    
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
static std::vector<grid_cell_t *> rank_ghost_grid_cells;
static std::vector<double> bottom_left_point;
static std::vector<double> top_right_point;

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

// Helper function for returning the rank of that a particle cell belongs.
int get_rank_from_particle_position(particle_t* p) {
    return grid_cells[get_block_index(p)]->rank;
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

    /* printf("Grid_dimension: %d\n", grid_dimension);
    printf("grid_cell_length: %f\n", grid_cell_length);
    printf("total_grid_count: %d\n", grid_dimension);
    fflush(stdout); */

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
    //printf("Size: %f\n", size);
    //fflush(stdout);

    // Assign grid_cell_t for each rank
    int starting_row = 0;
    for (int k = 0; k < num_procs; k += 1) {
        int rows_assigned_to_rank = grid_dimension / num_procs; // int division
        if (grid_dimension % num_procs > k ) { //n rows leftover, assign to first n ranks
            rows_assigned_to_rank += 1;
        }
        //printf("Starting_row: %d, k: %d, Rank: %d\n",starting_row, k, rank);
        //fflush(stdout);

        for (int i = 0; i < rows_assigned_to_rank; i++) {
            for (int j = 0; j < grid_dimension; j++) {
                // Compute current grid index and get the corresponding grid_cell_t
                int curr_grid_index = (i + starting_row) * grid_dimension + j;
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
        starting_row += rows_assigned_to_rank; // incrementing for next rank's iteration
        
    }

    //TODO: consider if we need the ghost_particles attribute
    // Assign ghost particles AND update is_ghost_cell_to
    for (int i = 0; i < grid_cells.size(); i += 1) {
        grid_cell_t* curr_grid_cell = grid_cells[i];
        // Up
        if (curr_grid_cell->neighbor_up != NULL && curr_grid_cell->neighbor_up->rank != curr_grid_cell->rank) {

            curr_grid_cell->is_ghost_cell_to.insert(curr_grid_cell->neighbor_up->rank);
            curr_grid_cell->neighbor_up->is_ghost_cell_to.insert(curr_grid_cell->rank);
            rank_ghost_grid_cells.push_back(curr_grid_cell->neighbor_up);
        }
        // Down
        if (curr_grid_cell->neighbor_down != NULL && curr_grid_cell->neighbor_down->rank != curr_grid_cell->rank) {

            curr_grid_cell->is_ghost_cell_to.insert(curr_grid_cell->neighbor_down->rank);
            curr_grid_cell->neighbor_down->is_ghost_cell_to.insert(curr_grid_cell->rank);
            rank_ghost_grid_cells.push_back(curr_grid_cell->neighbor_down);
        }
        // Left
        if (curr_grid_cell->neighbor_left != NULL && curr_grid_cell->neighbor_left->rank != curr_grid_cell->rank) {

            curr_grid_cell->is_ghost_cell_to.insert(curr_grid_cell->neighbor_left->rank);
            curr_grid_cell->neighbor_left->is_ghost_cell_to.insert(curr_grid_cell->rank);
            rank_ghost_grid_cells.push_back(curr_grid_cell->neighbor_left);
        }
        // Right
        if (curr_grid_cell->neighbor_right != NULL && curr_grid_cell->neighbor_right->rank != curr_grid_cell->rank) {

            curr_grid_cell->is_ghost_cell_to.insert(curr_grid_cell->neighbor_right->rank);  
            curr_grid_cell->neighbor_right->is_ghost_cell_to.insert(curr_grid_cell->rank);
            rank_ghost_grid_cells.push_back(curr_grid_cell->neighbor_right);
        }
        // Up-Left
        if (curr_grid_cell->neighbor_left_up != NULL && curr_grid_cell->neighbor_left_up->rank != curr_grid_cell->rank) {

            curr_grid_cell->is_ghost_cell_to.insert(curr_grid_cell->neighbor_left_up->rank);
            curr_grid_cell->neighbor_left_up->is_ghost_cell_to.insert(curr_grid_cell->rank);
            rank_ghost_grid_cells.push_back(curr_grid_cell->neighbor_left_up);
        }
        // Up-Right
        if (curr_grid_cell->neighbor_right_up != NULL && curr_grid_cell->neighbor_right_up->rank != curr_grid_cell->rank) {
            curr_grid_cell->is_ghost_cell_to.insert(curr_grid_cell->neighbor_right_up->rank);
            curr_grid_cell->neighbor_right_up->is_ghost_cell_to.insert(curr_grid_cell->rank);
            rank_ghost_grid_cells.push_back(curr_grid_cell->neighbor_right_up);
        }
        // Down-Left
        if (curr_grid_cell->neighbor_left_down != NULL && curr_grid_cell->neighbor_left_down->rank != curr_grid_cell->rank) {
            curr_grid_cell->is_ghost_cell_to.insert(curr_grid_cell->neighbor_left_down->rank);
            curr_grid_cell->neighbor_left_down->is_ghost_cell_to.insert(curr_grid_cell->rank);
            rank_ghost_grid_cells.push_back(curr_grid_cell->neighbor_left_down);
        }
        // Down-Right
        if (curr_grid_cell->neighbor_right_down != NULL && curr_grid_cell->neighbor_right_down->rank != curr_grid_cell->rank) {
            //Updates is_ghost_cell_to
            curr_grid_cell->is_ghost_cell_to.insert(curr_grid_cell->neighbor_right_down->rank);
            curr_grid_cell->neighbor_right_down->is_ghost_cell_to.insert(curr_grid_cell->rank);
            rank_ghost_grid_cells.push_back(curr_grid_cell->neighbor_right_down);
        }
    }

    /* for (int i = 0; i < grid_cells.size(); i += 1) {
        grid_cell_t* curr_grid_cell = grid_cells[i];
        printf("Rank: %d, cell_id: %d, is_ghost_cell_to count: %d\n", 
        rank, curr_grid_cell->index, curr_grid_cell->is_ghost_cell_to.size());
        if (curr_grid_cell->is_ghost_cell_to.size() > 0) {
            printf("Ghost cell of rank: %d\n", *curr_grid_cell->is_ghost_cell_to.begin());
        }
    }
    fflush(stdout); */
    
    // More logging to confirm end of fxn was reached
    //printf("made it to end of init_simulation\n");
    //fflush(stdout);
}

// Return a vector of neighbor grid cells. Self-exclusive.
std::vector<grid_cell_t*> get_neighbor_cells(grid_cell_t* grid) {
    std::vector<grid_cell_t*> neighbor_vector;
    // left-up
    if (grid->neighbor_left_up) {
        neighbor_vector.push_back(grid->neighbor_left_up);
    }
    // up
    if (grid->neighbor_up) {
        neighbor_vector.push_back(grid->neighbor_up);
    }
    // right-up
    if (grid->neighbor_right_up) {
        neighbor_vector.push_back(grid->neighbor_right_up);
    }
    // left
    if (grid->neighbor_left) {
        neighbor_vector.push_back(grid->neighbor_left);
    }
    // right
    if (grid->neighbor_right) {
        neighbor_vector.push_back(grid->neighbor_right);
    }
    // left-down
    if (grid->neighbor_left_down) {
        neighbor_vector.push_back(grid->neighbor_left_down);
    }
    // down
    if (grid->neighbor_down) {
        neighbor_vector.push_back(grid->neighbor_down);
    }
    // right-down
    if (grid->neighbor_right_down) {
        neighbor_vector.push_back(grid->neighbor_right_down);
    }
    return neighbor_vector;
}


void simulate_one_step(particle_t* parts, int num_parts, double size, int rank, int num_procs, int step) {
    // Write this function

    // For each grid cell object calculate forces between all particles
        // For particles within current grid cell
            // Intercel force 
                // Apply force
            // Intracel force
                // For neighbor particle
                    // Apply force
    
    for (grid_cell_t* g: rank_grid_cells) {
        for (particle_t* p: g->particles) {
            
            p->ax = p->ay = 0;
            //Forces within the cell
            for (particle_t* n: g->particles) {
                apply_force(*p, *n);
                
            }
    

            //Forces from neighbor cells
            std::vector<grid_cell_t*> neighbor_grid_cells = get_neighbor_cells(g);
            
            for (grid_cell_t* neighbor_g: neighbor_grid_cells){
                for (particle_t* n: neighbor_g->particles) {
                    apply_force(*p, *n);
                    
                }
            }
        }
    }

    //printf("Done applying forces, rank: %d\n", rank);
    //fflush(stdout);

    // Save particles to send (by p.ID), for each neighbor rank
    std::vector<std::unordered_set<u_int64_t>> particle_ids_to_send;
    particle_ids_to_send.resize(num_procs);

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
            //       The 3 conditions can occur simultaneously 
  	for (grid_cell_t* g: rank_grid_cells) {
        // TODO: simple optimization, maybe have a list of such grid_cells
        if (g->is_ghost_cell_to.size() > 0) {
            // Case 2: all particles that started in this grid needs to be communicated
            for (particle_t* p: g->particles){
                for (int neighbor_rank: g->is_ghost_cell_to) {
                    particle_ids_to_send[neighbor_rank].insert(p->id);
                }
            }
        } 
    }

    //printf("Done checking ghost_grids, rank: %d\n", rank);
    //fflush(stdout);

    // A vector of ids of particles owned by this rank. Later used for grid updates.
    std::vector<int> rank_part_ids_before_move;

    for (grid_cell_t* g: rank_grid_cells) {
        for (particle_t* p: g->particles) {
            rank_part_ids_before_move.push_back(p->id);
            move(*p, size);

            int new_rank = get_rank_from_particle_position(p);
            if (new_rank != rank) { // Case 1
                particle_ids_to_send[new_rank].insert(p->id);
            }
            grid_cell_t* grid_cell_after_move = grid_cells[get_block_index(p)];
            if (grid_cell_after_move->is_ghost_cell_to.size() > 0) { // Case 3
                for (int neighbor_rank: grid_cell_after_move->is_ghost_cell_to) {
                    particle_ids_to_send[neighbor_rank].insert(p->id);
                }
            }
        }
    }
    
    // Update grids for particles owned by the rank pre-move
    for (grid_cell_t* g: grid_cells) {
        g->particles.clear();
    }
    for (int part_id: rank_part_ids_before_move) {
        particle_t* p = parts + (part_id - 1);
        grid_cells[get_block_index(p)]->particles.push_back(p);
        // TODO: Update if:
        // 1. Still in chunk owned by this rank
        // 2. If in ghost zone of this rank.
    }


  	// Communication section
  		// For each neighbor rank
  			// Put particles meant for that rank into buffer, send to that rank
  			// MPI_Probe the amount of particles to receive, allocate buffer, receive particles
  			// update local particles in PARTS

    // First exchange the amount of particles to send/receive
    int* particle_counts_to_send;
    int* particle_counts_to_receive;

    MPI_Alloc_mem(num_procs * sizeof(int), MPI_INFO_NULL, &particle_counts_to_send);
    MPI_Alloc_mem(num_procs * sizeof(int), MPI_INFO_NULL, &particle_counts_to_receive);
    for (int target_rank = 0; target_rank < num_procs; target_rank++) {
        if (target_rank == rank) { // No need to send to self
            particle_counts_to_send[target_rank] = 0;
        } else {
            particle_counts_to_send[target_rank] = particle_ids_to_send[target_rank].size();
            //printf("Rank: %d, target_rank: %d, size: %d\n", 
            //rank, target_rank, particle_counts_to_send[target_rank]);
            //fflush(stdout);
        }
    }
    //printf("Counts_to_send[0]: %d, Counts_to_send[1]: %d, Rank: %d\n",
    //    particle_counts_to_send[0], particle_counts_to_send[1], rank);
    //printf("Before_alltoall, rank: %d\n", rank);
    //fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);
    // Need to barrier here? - is there a better way?

    // Use all_to_all, getting the amount of particles to receive from each rank
    MPI_Alltoall(particle_counts_to_send, 1, MPI_INT, \
        particle_counts_to_receive, 1, MPI_INT, MPI_COMM_WORLD);

    //printf("After_alltoall, rank: %d\n", rank);
    //fflush(stdout);
    // Prepare send and receive buffers
    std::vector<std::vector<particle_t>> send_buffers; // keys = rank to send to
    std::vector<std::vector<particle_t>> recv_buffers; // keys = rank to receive from
    send_buffers.resize(num_procs);
    recv_buffers.resize(num_procs);
    // Fill send buffers, resize receive buffers
    for (int r = 0; r < num_procs; r++) {
        for (int part_id: particle_ids_to_send[r]) {
            send_buffers[r].push_back(parts[part_id - 1]);
        }
        recv_buffers[r].resize(particle_counts_to_receive[r]);
    }

    // Async sends and receives
    MPI_Request send_requests[num_procs];
    MPI_Request recv_requests[num_procs];
    for (int r = 0; r < num_procs; r++) {
        MPI_Isend(&send_buffers[r][0], send_buffers[r].size(), PARTICLE, \
            r, 0, MPI_COMM_WORLD, &send_requests[r]);
        MPI_Irecv(&recv_buffers[r][0], recv_buffers[r].size(), PARTICLE, \
            r, 0, MPI_COMM_WORLD, &recv_requests[r]);
    }

    MPI_Status send_statuses[num_procs];
    MPI_Status recv_statuses[num_procs];
    // Wait on all the requests
    MPI_Waitall(num_procs, send_requests, send_statuses);
    MPI_Waitall(num_procs, recv_requests, recv_statuses);
    
    // For each received particle:
    //  Update received particles info to local copy in PARTS
    //  Place into rank-local corresponding grid_cells 
    for (int r = 0; r < num_procs; r++) {
        for (particle_t p: recv_buffers[r]) {
            parts[p.id - 1] = p;
            int grid_index = get_block_index(&p);
            grid_cells[grid_index]->particles.push_back(&parts[p.id - 1]);
        }
    }

    MPI_Free_mem(particle_counts_to_send);
    MPI_Free_mem(particle_counts_to_receive);

}

void gather_for_save(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    // Write this function such that at the end of it, the master (rank == 0)
    // processor has an in-order view of all particles. That is, the array
    // parts is complete and sorted by particle id.

    // Log how many particles were passerd in for debugging
    // for (int i = 0; i < num_parts; i += 1) {
    //     printf("id: %d\n", parts[i].id);
    // }
    //printf("made it to beggining of gather_for_save\n");
    //fflush(stdout);

    // Vectors for gathering particles from processors
    std::vector<particle_t> sending_parts;  // Size depends number of particles owned by each processor
    std::vector<particle_t> receiving_parts(num_parts); // Size

    // Add particles from processor to be sent for gathering
    for (int i = 0; i < rank_grid_cells.size(); i += 1) {
        std::vector<particle_t*> curr_grid_cell_parts = rank_grid_cells[i]->particles;
        for (int j = 0; j < curr_grid_cell_parts.size(); j += 1) {
            sending_parts.push_back(*curr_grid_cell_parts[j]);
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
    //printf("rank %d before MPI_Gather\n", rank);
    //fflush(stdout);
    int sending_count = sending_parts.size();
    MPI_Gather(&sending_count, 1, MPI_INT, receiving_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Build up offsets array
    if (rank == 0) {
        offsets[0] = 0;
        for (int i = 1; i < num_procs; i += 1) {
            offsets[i] = offsets[i - 1] + receiving_counts[i - 1];
        }
        //printf("rank: %d, receiving_counts: %d, %d\n", rank, receiving_counts[0], receiving_counts[1]);
        //printf("rank: %d, offsets: %d, %d\n", rank, offsets[0], offsets[1]);
    }

    // Variable gather of particles from processors
    //printf("rank %d before MPI_Gatherv\n", rank);
    //printf("rank: %d, size: %d\n", rank, sending_parts.size());   
    //fflush(stdout);
    MPI_Gatherv(&sending_parts[0], sending_parts.size(), PARTICLE, &receiving_parts[0], receiving_counts, offsets, PARTICLE, 0, MPI_COMM_WORLD);

    //printf("rank %d after MPI_Gatherv\n", rank);
    //fflush(stdout);
    // Create in-order view of all particles, sorted by particle id
    if (rank == 0) {
        for (int i = 0; i < num_parts; i += 1) {
            particle_t curr_part = receiving_parts[i];
            parts[curr_part.id - 1] = curr_part;
        }
    }

    //printf("made it to end of gather_for_save\n");
    //fflush(stdout);
}