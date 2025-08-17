#include <iostream>
#include <vector>
#include <string>
#include <map>

#include <chrono>
#include <limits>



#include <ginkgo/ginkgo.hpp>

#include "geometry.hpp"
#include "IO_VTK.hpp"
#include "geo_math.hpp"

#include "MPI_utils.hpp"

#include <omp.h>
#include <mpi.h>


int main(int argc, char* argv[]) {


    MPI_Init(NULL, NULL);

    MPI_Status status;



    int my_rank, world_size;
    int ierror;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);




    create_MPIFace_type();
    create_MPIVertex_type();
    int dynamic_omp_nthreads{0};

    int total_n_vertices{0}, total_n_faces{0};
    int local_n_vertices{0}, local_n_faces{0};
    int faces_base{0}, faces_remainder{0};
    int faces_chunk_size{0};
    int max_chunk_faces{0}, max_chunk_vertices{0};
    int n_ownrank_repeating_vertices{0};
    int distribute_faces_offset{0};

    int curr_block_offset{0};
    // height constant per rank, only width block changes
    int curr_block_height{0};
    int curr_block_width{0};

    int target_ring_rank{0};
    int source_ring_rank{0};

    double max_construction_time{0.0};

    float cent[3];
    
    float* G_arr;

    std::vector<Vertex> all_vertices;
    std::vector<Face> all_faces;
    std::vector<Vertex> local_vertices;
    std::vector<Vertex> repeating_vertices_send_buffer;
    std::vector<Face> local_faces;

    std::map<int, int> global_id_to_local;


    std::string mesh_file = (argc > 1 ? argv[1] : "sphere");
    
    dynamic_omp_nthreads = omp_get_thread_num();

    if (my_rank == 0) {

        read_vtk(mesh_file, all_vertices, all_faces, total_n_vertices, total_n_faces);
//        std::cout << "Prod. 1 finished reading the geometry!" << std::endl 
//                  << "Mesh: " << mesh_file << ", " << total_n_vertices << 
 //                    " nodes, " << total_n_faces << " faces" << std::endl;

        G_arr = (float *) calloc(total_n_faces * total_n_faces, sizeof(float));
    }
      


    MPI_Bcast(&total_n_faces, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&total_n_vertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);


    // The way people do it in mpi

    faces_base = total_n_faces / world_size;
    faces_remainder = total_n_faces % world_size;

    local_n_faces = (my_rank < faces_remainder) ? faces_base + 1 : faces_base;

    max_chunk_faces = faces_base + 1;

    local_faces.resize(max_chunk_faces);


    if (my_rank == 0) {

        int n_repeating_vertices_per_rank{0};

        std::copy(all_faces.begin(), all_faces.begin() + local_n_faces, local_faces.begin());


        n_ownrank_repeating_vertices = local_n_faces * 3;
        local_vertices.reserve(n_ownrank_repeating_vertices);

        // TODO : if this is too slow, move later
        for (int j = 0; j < local_n_faces; j++) {
                local_vertices.push_back(Vertex(all_vertices[all_faces[j].v1]));
                local_vertices.push_back(Vertex(all_vertices[all_faces[j].v2]));
                local_vertices.push_back(Vertex(all_vertices[all_faces[j].v3]));
        }


        distribute_faces_offset = local_n_faces;

        for (int dest_rank = 1; dest_rank < world_size; dest_rank++) {

            faces_chunk_size = (dest_rank < faces_remainder) ? faces_base + 1 : faces_base;

            //forward the faces
            MPI_Send(all_faces.data() + distribute_faces_offset,
                     faces_chunk_size, MPI_FACE, dest_rank, 44, MPI_COMM_WORLD);

            n_repeating_vertices_per_rank = faces_chunk_size * 3;
            repeating_vertices_send_buffer.clear();
            repeating_vertices_send_buffer.reserve(n_repeating_vertices_per_rank);

            for (int j = distribute_faces_offset; j < distribute_faces_offset + faces_chunk_size; j++) {

                repeating_vertices_send_buffer.push_back(Vertex(all_vertices[all_faces[j].v1]));
                repeating_vertices_send_buffer.push_back(Vertex(all_vertices[all_faces[j].v2]));
                repeating_vertices_send_buffer.push_back(Vertex(all_vertices[all_faces[j].v3]));

            }


            MPI_Send(repeating_vertices_send_buffer.data(),
                n_repeating_vertices_per_rank, MPI_VERTEX, dest_rank, 43, MPI_COMM_WORLD);


            distribute_faces_offset += faces_chunk_size;

        }

        repeating_vertices_send_buffer.clear();

    }

    
    else {






        MPI_Recv(local_faces.data(), local_n_faces, MPI_FACE, 0, 44, MPI_COMM_WORLD, &status);


        n_ownrank_repeating_vertices = local_n_faces * 3;
        local_vertices.resize(n_ownrank_repeating_vertices);

        MPI_Recv(local_vertices.data(), n_ownrank_repeating_vertices, 
            MPI_VERTEX, 0, 43, MPI_COMM_WORLD, &status);



    }

    

    std::sort(local_vertices.begin(), local_vertices.end());
    auto last = std::unique(local_vertices.begin(), local_vertices.end());
    local_vertices.erase(last, local_vertices.end());

    local_n_vertices = local_vertices.size();
    

    // to share 

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(&local_n_vertices, &max_chunk_vertices, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    local_vertices.resize(max_chunk_vertices);



    // map vertex ids to local ids
    for (int i = 0; i < local_n_vertices; i++) {
        global_id_to_local[local_vertices[i].id] = i;
    }

    for (int i = 0; i < local_n_faces; i++) {
        local_faces[i].v1 = global_id_to_local[local_faces[i].v1];
        local_faces[i].v2 = global_id_to_local[local_faces[i].v2];
        local_faces[i].v3 = global_id_to_local[local_faces[i].v3];
    }

    global_id_to_local.clear();

    //switching to C style loops for convenience when loading to GPU
    // of size max_chunk faces to also be able later receive buffers of other ranks
    float (*centroids)[3] = (float (*)[3])calloc(max_chunk_faces, sizeof(float[3]));


    for (int i = 0; i < local_n_faces; ++i) {
        centroids[i][0] = (local_vertices[local_faces[i].v1].x + local_vertices[local_faces[i].v2].x + local_vertices[local_faces[i].v3].x) / 3;
        centroids[i][1] = (local_vertices[local_faces[i].v1].y + local_vertices[local_faces[i].v2].y + local_vertices[local_faces[i].v3].y) / 3;
        centroids[i][2] = (local_vertices[local_faces[i].v1].z + local_vertices[local_faces[i].v2].z + local_vertices[local_faces[i].v3].z) / 3;
    }

    float* G_arr_block = (float *) calloc(total_n_faces * local_n_faces, sizeof(float));



    MPI_Exscan(&local_n_faces, &curr_block_offset,
               1, MPI_INT, MPI_SUM,
               MPI_COMM_WORLD);
                                                                     

    curr_block_height = local_n_faces;
    curr_block_width = local_n_faces;

    double start_time{0.0}, end_time{0.0}, assem_time{0.0};
    start_time = MPI_Wtime();
    for (int rotation_count = 0; rotation_count < world_size; rotation_count++) {
        
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < curr_block_height; i++) {
            for (int j = 0; j < curr_block_width; j++) {

                        G_arr_block[total_n_faces * i + curr_block_offset + j] =
                            ( (i == j) && (rotation_count == 0) )
                                ? regularized_integral(
                                      local_vertices[local_faces[i].v1],
                                      local_vertices[local_faces[i].v2],
                                      local_vertices[local_faces[i].v3])
                                : gauss_integral(
                                      local_vertices[local_faces[j].v1],
                                      local_vertices[local_faces[j].v2],
                                      local_vertices[local_faces[j].v3],
                                      centroids[i]);

            }
        }


        int target_ring_rank = (my_rank - 1 + world_size) % world_size; 
        int source_ring_rank = (my_rank + 1) % world_size; 


        curr_block_offset = (curr_block_offset + curr_block_width) % total_n_faces;

        MPI_Sendrecv_replace(&curr_block_width, 1, MPI_INT, target_ring_rank, rotation_count,  
            source_ring_rank, rotation_count, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv_replace(local_faces.data(), max_chunk_faces, MPI_FACE, target_ring_rank,
           100 + rotation_count, source_ring_rank,   
           100 + rotation_count,
            MPI_COMM_WORLD, &status
        );

        MPI_Sendrecv_replace(local_vertices.data(), max_chunk_vertices, MPI_VERTEX, target_ring_rank,
           200 + rotation_count,
           source_ring_rank,   200 + rotation_count,
            MPI_COMM_WORLD, &status
        );


    }

    end_time = MPI_Wtime();
    assem_time = end_time - start_time;
    MPI_Barrier(MPI_COMM_WORLD);

    

    MPI_Reduce(&assem_time, &max_construction_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);




    // Gather send counts (in floats) from each rank
    int block_volume = curr_block_height * total_n_faces; // number of floats per rank
    int offset_volume = curr_block_offset * total_n_faces; // offset per number of floats

    // time to go to big data, not very pleasent
    std::vector<int> recvcounts(world_size);
    std::vector<int> offsets_in_rows(world_size);

    MPI_Gather(&block_volume, 1, MPI_INT,
               recvcounts.data(), 1, MPI_INT,
               0, MPI_COMM_WORLD);

    // Gather curr_block_offset from each rank (in rows)
    MPI_Gather(&offset_volume, 1, MPI_INT,
               offsets_in_rows.data(), 1, MPI_INT,
               0, MPI_COMM_WORLD);

    // Build displacements (in floats) from offsets_in_rows

    // Now gather the actual data into G_arr on root
    MPI_Gatherv(G_arr_block, block_volume, MPI_FLOAT,
                my_rank == 0 ? G_arr : nullptr,
                recvcounts.data(), offsets_in_rows.data(), MPI_FLOAT,
                0, MPI_COMM_WORLD);



    local_vertices.clear();
    local_faces.clear();
    free(G_arr_block);

    if (my_rank == 0) {

    std::cout <<  max_construction_time;

    // proceed with the rest of the computation: serialized
    float* phi_arr = (float *) calloc(total_n_faces, sizeof(float));
    float* q_arr = (float *) calloc(total_n_faces, sizeof(float));


    // Boundary Conditions
    for (auto& v :all_vertices) { v.potential = v.x * v.x + v.y * v.y;} 
    
    for (size_t i = 0; i <total_n_faces; ++i) {
        phi_arr[i] = (all_vertices[all_faces[i].v1].potential + all_vertices[all_faces[i].v2].potential + all_vertices[all_faces[i].v3].potential) / 3.0;
    }



  //  std::cout << gko::version_info::get() << std::endl;

    omp_set_num_threads(32); 
    const auto exec = gko::OmpExecutor::create();


    //system matrix
    std::shared_ptr<gko::matrix::Dense<float>> G = 
            gko::matrix::Dense<float>::create(                             
                exec, gko::dim<2>{total_n_faces,total_n_faces},      
                gko::array<float>::view(exec,total_n_faces *total_n_faces, G_arr)
                ,total_n_faces); 

    //RHS, load vector
    std::shared_ptr<gko::matrix::Dense<float>> phi = 
            gko::matrix::Dense<float>::create(                             
                exec, gko::dim<2>{total_n_faces, 1},      
                gko::array<float>::view(exec,total_n_faces, phi_arr)
                , 1); 

    // to store solution
    std::shared_ptr<gko::matrix::Dense<float>> q = 
            gko::matrix::Dense<float>::create(                             
                exec, gko::dim<2>{total_n_faces, 1},      
                gko::array<float>::view(exec,total_n_faces, q_arr)
                , 1);

    auto gmres_gen = gko::solver::Gmres<float>::build()
        .with_criteria(
            gko::stop::Iteration::build().with_max_iters(200u).on(exec),
            gko::stop::ResidualNorm<float>::build()
                .with_reduction_factor(1e-6).on(exec))
        .with_preconditioner(
            gko::preconditioner::Jacobi<float, int>::build()
                .with_max_block_size(1u)
                .on(exec))
        .on(exec);



    auto gmres_solver = gmres_gen->generate(G);



    // solver call

    double start_solveTime{0.0}, end_solveTime{0.0};
    start_solveTime = MPI_Wtime();
    gmres_solver->apply(phi, q);
    end_solveTime = MPI_Wtime();

   // std::cout << "Time to solve: " << end_solveTime - start_solveTime << "s" << std::endl;

    free(G_arr);
    free(phi_arr);

    std::vector vertex_densities(total_n_vertices, 0.0); 
    std::vector ring_areas(total_n_vertices, 0.0);

    for (int i = 0; i < total_n_faces; ++i) {
        auto& f = all_faces[i];
        double Ai = face_area(all_vertices[f.v1], all_vertices[f.v2], all_vertices[f.v3]);
        
        vertex_densities[all_vertices[f.v1].id] += q_arr[i] * Ai;
        ring_areas[all_vertices[f.v1].id] += Ai;

        vertex_densities[all_vertices[f.v2].id] += q_arr[i] * Ai;
        ring_areas[all_vertices[f.v2].id] += Ai;

        vertex_densities[all_vertices[f.v3].id] += q_arr[i] * Ai;
        ring_areas[all_vertices[f.v3].id] += Ai;
    }

    for (int i = 0; i < total_n_vertices; ++i) {
        if (ring_areas[i] > 1e-12 /* diving by very small areas */) { 
            vertex_densities[i] /= ring_areas[i];
        }
    }
    ring_areas.clear();
    free(q_arr);



    std::vector<std::vector<int>> adjacency(total_n_vertices);
    for (const auto& face : all_faces) {
        int v1_id = all_vertices[face.v1].id;
        int v2_id = all_vertices[face.v2].id;
        int v3_id = all_vertices[face.v3].id;
    
	// Add edges to adjacency list, avoiding duplicates
        adjacency[v1_id].push_back(v2_id); adjacency[v1_id].push_back(v3_id);
        adjacency[v2_id].push_back(v1_id); adjacency[v2_id].push_back(v3_id);
        adjacency[v3_id].push_back(v1_id); adjacency[v3_id].push_back(v2_id);
    }
    for(auto& neighbors : adjacency) {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }



    const int smoothing_iterations = 1; 
    std::vector<double> smoothed_densities = vertex_densities; 

    for (int iter = 0; iter < smoothing_iterations; ++iter) {
        for (int i = 0; i < total_n_vertices; ++i) {
            if (adjacency[i].empty()) continue;

            double neighbor_sum = 0.0;
            for (int neighbor_id : adjacency[i]) {
                neighbor_sum += vertex_densities[neighbor_id];
            }
            smoothed_densities[i] = neighbor_sum / adjacency[i].size();
        }
        vertex_densities = smoothed_densities; 
    }


    for (auto& v : all_vertices) {
        v.density = vertex_densities[v.id];
    }

   // std::cout << "Smoothing complete." << std::endl;



    write_vtu(mesh_file, all_vertices, all_faces, total_n_vertices, total_n_faces);
    }
    MPI_Finalize();

    return 0;

}

