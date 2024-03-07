#include "common.h"
#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <set>
#include <algorithm>

/*
* Static global variables
*/
struct grid_cell_t {
    int rank;
    int index;
    std::vector<particle_t*> particles;
    std::set<int> is_ghost_cell_to; // stores rank indices
    std::vector<grid_cell_t*> neighbor_grids;

    // Constructor to initialize member vectors
    grid_cell_t() 
        : rank(0), 
          index(0),
          particles(),  // Initialize particle vector (empty)
          is_ghost_cell_to(), // Initialize is_ghost_cell_to vector (empty)
          neighbor_grids() {} 
};

static int grid_dimension;
static double grid_cell_length;
static std::vector<grid_cell_t> grid_cells; // "source truth"
static std::vector<grid_cell_t*> rank_grid_cells;
static std::unordered_set<grid_cell_t *> rank_ghost_grid_cells;
static std::vector<grid_cell_t *> rank_grid_cells_that_are_ghost;
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
    return grid_cells[get_block_index(p)].rank;
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

    printf("Grid_dimension: %d\n", grid_dimension);
    printf("grid_cell_length: %f\n", grid_cell_length);
    printf("total_grid_count: %d\n", grid_dimension);
    fflush(stdout);

    // Create grid of grid_cell_t
    grid_cells.resize(total_grid_count);
    /* for (int i = 0; i < total_grid_count; ++i) {
        // Create a new grid_cell_t instance and assign its pointer to the vector
        grid_cells[i] = new grid_cell_t();
    } */

    // Iterate over all pariticles and place them into their corresponding grid_cell_t
    for (int i = 0; i < num_parts; i++) {
        particle_t* particle = parts + i;
        int block_index = get_block_index(particle);
        grid_cells[block_index].particles.push_back(particle);
    }

    // Assign grid_cell_t for each rank
    // TODO: this needs to change for different domain division methods
    int starting_row = 0;
    for (int k = 0; k < num_procs; k += 1) {
        int rows_assigned_to_rank = grid_dimension / num_procs; // int division
        if (grid_dimension % num_procs > k ) { //n rows leftover, assign to first n ranks
            rows_assigned_to_rank += 1;
        }

        for (int i = 0; i < rows_assigned_to_rank; i++) {
            for (int j = 0; j < grid_dimension; j++) {
                // Compute current grid index and get the corresponding grid_cell_t
                int curr_grid_index = (i + starting_row) * grid_dimension + j;
                grid_cell_t curr_grid_cell = grid_cells[curr_grid_index];

                // Assign rank and index to current grid_cell
                curr_grid_cell.rank = k;
                curr_grid_cell.index = curr_grid_index;

                // Add grid_cell to this rank's vector of grid_cell_t
                if (rank == k) {
                    rank_grid_cells.push_back(&grid_cells[curr_grid_index]);
                }
            }
        }
        starting_row += rows_assigned_to_rank; // incrementing for next rank's iteration 
    }

    //"linking" grid_cells
    for (int g = 0; g < total_grid_count; g++) {
        int col_index = g % grid_dimension;
        int row_index = g / grid_dimension;
        int row_min = fmax(0, row_index - 1);
        int row_max = fmin(grid_dimension - 1, row_index + 1);
        int col_min = fmax(0, col_index - 1);
        int col_max = fmin(grid_dimension - 1, col_index + 1);

        if (col_index < grid_dimension - 1) {
            for (int row = row_min; row <= row_max; row++) {
                int neighbor_grid_index = col_index + 1 + row * grid_dimension;
                
                grid_cells[g].neighbor_grids.push_back(&grid_cells[neighbor_grid_index]);
                grid_cells[neighbor_grid_index].neighbor_grids.push_back(&grid_cells[g]);
            }
        }
        if (row_index > 0) {
            int neighbor_grid_index = g - grid_dimension;
            grid_cells[g].neighbor_grids.push_back(&grid_cells[neighbor_grid_index]);
            grid_cells[neighbor_grid_index].neighbor_grids.push_back(&grid_cells[g]);
        }
    }

    //Updating is_ghost_to
    for (int g = 0; g < total_grid_count; g++) {
        grid_cell_t curr_grid_cell = grid_cells[g];
        for (grid_cell_t* neighbor_grid_cell: curr_grid_cell.neighbor_grids) {
            if (curr_grid_cell.rank != neighbor_grid_cell->rank) {
                curr_grid_cell.is_ghost_cell_to.insert(neighbor_grid_cell->rank);

                //TODO: one way is ok if iterating through all grids, since neighborhood is always mutual
                //neighbor_g->is_ghost_cell_to.insert(curr_grid_cell->rank);
                rank_ghost_grid_cells.insert(neighbor_grid_cell);  
            }
        }
    }

    //Update rank_grid_cells_that_are_ghost (aka cells owned by this rank, and is a ghost cell to some other rank)
    for (grid_cell_t* g: rank_grid_cells) {
        if (g->is_ghost_cell_to.size() > 0) {
            rank_grid_cells_that_are_ghost.push_back(g);
        }
    }
    
    // More logging to confirm end of fxn was reached
    //printf("made it to end of init_simulation\n");
    //fflush(stdout);
}

// Return a vector of neighbor grid cells. Self-exclusive.
std::vector<grid_cell_t*> get_neighbor_cells(grid_cell_t* grid) {
    std::vector<grid_cell_t*> neighbor_vector;
    
    return grid->neighbor_grids;
}

//TODO: change signature back
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
            for (grid_cell_t* neighbor_g: g->neighbor_grids){
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


    for (grid_cell_t* g: rank_grid_cells_that_are_ghost) {
        // Case 2: all particles that started in this grid needs to be communicated
        for (particle_t* p: g->particles){
            for (int neighbor_rank: g->is_ghost_cell_to) {
                particle_ids_to_send[neighbor_rank].insert(p->id);
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
            grid_cell_t grid_cell_after_move = grid_cells[get_block_index(p)];
            if (grid_cell_after_move.is_ghost_cell_to.size() > 0) { // Case 3
                for (int neighbor_rank: grid_cell_after_move.is_ghost_cell_to) {
                    particle_ids_to_send[neighbor_rank].insert(p->id);
                }
            }
        }
    }
    
    // Clear relevant grid_cells
    for (grid_cell_t* g: rank_grid_cells) {
        g->particles.clear();
    }
    for (grid_cell_t* g: rank_ghost_grid_cells) {
        g->particles.clear();
    }

    //printf("After clearing, rank: %d\n", rank);
    //fflush(stdout);

    // Update grids for particles owned by the rank pre-move
    for (int part_id: rank_part_ids_before_move) {
        particle_t* p = parts + (part_id - 1);

        // Place in our copy of grid if:
        // 1. Still in chunk owned by this rank
        // 2. If it is in ghost zone of this rank.
        grid_cell_t g = grid_cells[get_block_index(p)];
        if (g.rank == rank || rank_ghost_grid_cells.find(&g) != rank_ghost_grid_cells.end()) {
            grid_cells[get_block_index(p)].particles.push_back(p);
        }
    }

    //printf("Before communication section, rank: %d\n", rank);
    //fflush(stdout);

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
            particle_ids_to_send[target_rank].clear();
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

    //printf("After prepping buffers, rank: %d\n", rank);
    //fflush(stdout);

    // Async sends and receives

    //MPI_Request send_requests[num_procs];
    //MPI_Request recv_requests[num_procs];

    std::vector<MPI_Request> send_requests;
    std::vector<MPI_Request> recv_requests;

    int total_sends = 0;
    int total_recvs = 0;
    for (int r = 0; r < num_procs; r++) {
        if (particle_counts_to_send[r] > 0) {
            MPI_Request send_r;
            MPI_Isend(&send_buffers[r][0], send_buffers[r].size(), PARTICLE, \
                r, 0, MPI_COMM_WORLD, &send_r);
            send_requests.push_back(send_r);
            total_sends += 1;
        }
        if (particle_counts_to_receive[r] > 0) {
            MPI_Request recv_r;
            MPI_Irecv(&recv_buffers[r][0], recv_buffers[r].size(), PARTICLE, \
                r, 0, MPI_COMM_WORLD, &recv_r);
            recv_requests.push_back(recv_r);
            total_recvs += 1;
        }
    }

    //printf("After sends/recvs, rank: %d\n", rank);
    //fflush(stdout);

    MPI_Status send_statuses[total_sends];
    MPI_Status recv_statuses[total_recvs];
    
    // Wait on all the requests
    if (send_requests.size() > 0) {
        MPI_Waitall(total_sends, &send_requests[0], send_statuses);
    }
    if (recv_requests.size() > 0) {
        MPI_Waitall(total_recvs, &recv_requests[0], recv_statuses);
    }
    
    //printf("Before processing received parts, rank: %d\n", rank);
    //fflush(stdout);

    // For each received particle:
    //  Update received particles info to local copy in PARTS
    //  Place into rank-local corresponding grid_cells 
    for (int r = 0; r < num_procs; r++) {
        for (particle_t p: recv_buffers[r]) {
            parts[p.id - 1] = p;
            int grid_index = get_block_index(&p);
            grid_cells[grid_index].particles.push_back(&parts[p.id - 1]);
        }
    }
    //printf("After processing received parts, rank: %d\n", rank);
    //fflush(stdout);

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
        delete receiving_counts;
        delete offsets;
    }

    //printf("made it to end of gather_for_save\n");
    //fflush(stdout);
}