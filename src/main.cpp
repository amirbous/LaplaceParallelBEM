#include <iostream>
#include <vector>
#include <string>
#include <map>

#include <chrono>
#include <limits>


#include <ginkgo/ginkgo.hpp>

#include "../include/geometry.hpp"
#include "../include/IO_VTK.hpp"
#include "../include/geo_math.hpp"

#include "../include/MPI_utils.hpp"

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

    int total_n_vertices{0}, total_n_faces{0};
    int local_n_vertices{0}, local_n_faces{0};
    int faces_base{0}, faces_remainder{0};
    int faces_chunk_size{0};
    int max_chunk_faces{0}, max_chunk_vertices{0};
    int n_ownrank_repeating_vertices{0};
    int distribute_faces_offset{0};

    int curr_block_offset{0};

    std::vector<int> world_n_faces(world_size);
    std::vector<int> world_n_vertices(world_size);

    std::vector<Vertex> all_vertices;
    std::vector<Face> all_faces;
    std::vector<Vertex> local_vertices;
    std::vector<Vertex> repeating_vertices_send_buffer;
    std::vector<Face> local_faces;

    std::map<int, int> global_id_to_local;

    float cent[3];


    if (my_rank == 0) {

        

        std::string mesh_file = (argc > 1 ? argv[1] : "sphere");

        read_vtk(mesh_file, all_vertices, all_faces, total_n_vertices, total_n_faces);
        std::cout << "Prod. 1 finished reading the geometry!" << std::endl 
                  << "Mesh: " << mesh_file << ", " << total_n_vertices << 
                     " nodes, " << total_n_faces << " faces" << std::endl;
    }
      


    MPI_Bcast(&total_n_faces, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&total_n_vertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);


    // The way people do it in mpi

    faces_base = total_n_faces / world_size;
    faces_remainder = total_n_faces % world_size;

    local_n_faces = (my_rank < faces_remainder) ? faces_base + 1 : faces_base;

    max_chunk_faces = faces_base + 1;



    if (my_rank == 0) {

        int n_repeating_vertices_per_rank{0};
        world_n_faces[my_rank /* == 0 */] = local_n_faces;


        local_faces.assign(all_faces.begin(), all_faces.begin() + local_n_faces);



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
            world_n_faces[dest_rank] = faces_chunk_size;

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



        local_faces.resize(max_chunk_faces);


        MPI_Recv(local_faces.data(), local_n_faces, MPI_FACE, 0, 44, MPI_COMM_WORLD, &status);


        n_ownrank_repeating_vertices = local_n_faces * 3;
        local_vertices.resize(n_ownrank_repeating_vertices);

        MPI_Recv(local_vertices.data(), n_ownrank_repeating_vertices, 
            MPI_VERTEX, 0, 43, MPI_COMM_WORLD, &status);



    }

    
    MPI_Bcast(world_n_faces.data(), world_size, MPI_INT, 0, MPI_COMM_WORLD);


    std::sort(local_vertices.begin(), local_vertices.end());
    auto last = std::unique(local_vertices.begin(), local_vertices.end());
    local_vertices.erase(last, local_vertices.end());

    local_n_vertices = local_vertices.size();
    

    // to share 

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allgather(&local_n_vertices, 1, MPI_INT, world_n_vertices.data(), 1, MPI_INT, MPI_COMM_WORLD);

    MPI_Allreduce(&local_n_vertices, &max_chunk_vertices, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    local_vertices.resize(max_chunk_vertices);


    std::cout << "Proc. rank " << my_rank << " has " << local_n_faces << " faces and " << local_n_vertices << " vertices!" << std::endl;


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

    // precompute centroid which will be used later
    for (int i = 0; i < local_n_faces; ++i) {
        centroids[i][0] = (local_vertices[local_faces[i].v1].x + local_vertices[local_faces[i].v2].x + local_vertices[local_faces[i].v3].x) / 3;
        centroids[i][1] = (local_vertices[local_faces[i].v1].y + local_vertices[local_faces[i].v2].y + local_vertices[local_faces[i].v3].y) / 3;
        centroids[i][2] = (local_vertices[local_faces[i].v1].z + local_vertices[local_faces[i].v2].z + local_vertices[local_faces[i].v3].z) / 3;
    }

    float* G_arr = (float *) calloc(total_n_faces * local_n_faces, sizeof(float));


    // column offset for the diagonal
    curr_block_offset = std::accumulate(world_n_faces.begin(), world_n_faces.begin() + my_rank, 0);

    for (int i = 0; i < local_n_faces; i++) {
        G_arr[local_n_faces * i + curr_block_offset + i] = regularized_integral(local_vertices[local_faces[i].v1],
                                                                            local_vertices[local_faces[i].v2],
                                                                            local_vertices[local_faces[i].v3]);
    }



    
    for (int my_j = 0; my_j < world_size; my_j++) {
        if (my_j == my_rank) {
            std::cout << "rank " << my_j << std::endl;
            for (int i = 0; i < local_n_faces; i++) {

                std::cout << "(" << i + curr_block_offset <<  ", " << G_arr[local_n_faces * i + curr_block_offset + i] << ")" << std::endl;
            }
            std::cout << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);

    }








    MPI_Finalize();


    return 0;


}






    //MPI_Scatter(all_faces.data(), local_n_faces, MPI_FACE, local_faces.data(), local_n_faces,  
      //  MPI_FACE, 0, MPI_COMM_WORLD);




   // std::cout << std::endl;








/*


    ValueType* G_arr = (ValueType *) calloc(n_faces * n_faces, sizeof(ValueType));
    ValueType* phi_arr = (ValueType *) calloc(n_faces, sizeof(ValueType));
    ValueType* q_arr = (ValueType *) calloc(n_faces, sizeof(ValueType));


    // Boundary Conditions
    for (auto& v : vertices) { v.potential = 10;} 
    

    std::cout << "Started assembling the matrix" << std::endl;
    

   float (*centroids)[3] = (float (*)[3])malloc(n_faces * sizeof(float[3]));

    // Step 1: Compute all centroids in parallel
    for (size_t i = 0; i < n_faces; ++i) {
        centroids[i][0] = (vertices[faces[i].v1].x + vertices[faces[i].v2].x + vertices[faces[i].v3].x) / 3;
        centroids[i][1] = (vertices[faces[i].v1].y + vertices[faces[i].v2].y + vertices[faces[i].v3].y) / 3;
        centroids[i][2] = (vertices[faces[i].v1].z + vertices[faces[i].v2].z + vertices[faces[i].v3].z) / 3;
    }


    std::chrono::steady_clock::time_point begin_assembleMatrix = std::chrono::steady_clock::now(); // measure time for this section


    #pragma omp parallel for private(cent)
    for (size_t i = 0; i < n_faces; ++i) {
    // Compute centroid *once* for face i, private to each thread

    cent[0] = (vertices[faces[i].v1].x + vertices[faces[i].v2].x + vertices[faces[i].v3].x) / 3;
    cent[1] = (vertices[faces[i].v1].y + vertices[faces[i].v2].y + vertices[faces[i].v3].y) / 3;
    cent[2] = (vertices[faces[i].v1].z + vertices[faces[i].v2].z + vertices[faces[i].v3].z) / 3;

    for (size_t j = 0; j < n_faces; ++j) {


            G_arr[n_faces * i + j] = i == j ? 
                                        regularized_integral(vertices[faces[i].v1],
                                              vertices[faces[i].v2],
                                              vertices[faces[i].v3]) 
                                        :

                                        gauss_integral(vertices[faces[j].v1],
                                        vertices[faces[j].v2],
                                        vertices[faces[j].v3],
                                        cent);

        }
    }

    std::chrono::steady_clock::time_point end_assembleMatrix = std::chrono::steady_clock::now(); // end



    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < n_faces; ++i) {
        for (size_t j = 0; j < n_faces; ++j) {   
                G_arr[n_faces * i + j] = i == j ? 
                          regularized_integral(vertices[faces[i].v1],
                          vertices[faces[i].v2], vertices[faces[i].v3]) 
                                        :
                          gauss_integral(vertices[faces[j].v1], vertices[faces[j].v2],
                                         vertices[faces[j].v3], centroids[i]);
        }       
    }





    auto time_toAssemble = std::chrono::duration_cast<std::chrono::milliseconds>(end_assembleMatrix - begin_assembleMatrix).count();
    std::cout << "finished assembling the matrix in: " << time_toAssemble << "ms" << std::endl;

    // constructing load vector ==> density averaged over faces
    for (size_t i = 0; i < n_faces; ++i) {
        phi_arr[i] = (vertices[faces[i].v1].potential + vertices[faces[i].v2].potential + vertices[faces[i].v3].potential) / 3.0;
    }





    std::cout << gko::version_info::get() << std::endl;
    const auto exec = gko::OmpExecutor::create();

    //system matrix
    std::shared_ptr<gko::matrix::Dense> G = 
            gko::matrix::Dense::create(                             
                exec, gko::dim<2>{n_faces, n_faces},      
                gko::array::view(exec, n_faces * n_faces, G_arr)
                , n_faces); 

    //RHS, load vector
    std::shared_ptr<gko::matrix::Dense> phi = 
            gko::matrix::Dense::create(                             
                exec, gko::dim<2>{n_faces, 1},      
                gko::array::view(exec, n_faces, phi_arr)
                , 1); 

    // to store solution
    std::shared_ptr<gko::matrix::Dense> q = 
            gko::matrix::Dense::create(                             
                exec, gko::dim<2>{n_faces, 1},      
                gko::array::view(exec, n_faces, q_arr)
                , 1);

    auto gmres_gen = gko::solver::Gmres::build()
        .with_criteria(
            gko::stop::Iteration::build().with_max_iters(200u).on(exec),
            gko::stop::ResidualNorm::build()
                .with_reduction_factor(1e-6).on(exec))
        .with_preconditioner(
            gko::preconditioner::Jacobi<ValueType, int>::build()
                .with_max_block_size(1u)
                .on(exec))
        .on(exec);



    auto gmres_solver = gmres_gen->generate(G);


    std::chrono::steady_clock::time_point begin_solveMatrix = std::chrono::steady_clock::now(); // end
    
    // solver call
    gmres_solver->apply(phi, q);

    free(G_arr);
    free(phi_arr);

    std::chrono::steady_clock::time_point end_solveMatrix = std::chrono::steady_clock::now(); // end

    auto time_toSolve = std::chrono::duration_cast<std::chrono::milliseconds>(end_solveMatrix - begin_solveMatrix).count();
    
    std::cout << "Solved G*q = phi for charge distribution in: " << time_toSolve << "ms" << std::endl;


    // TODO: will be moved to a seperate function: considered not very important part of the code 
    // getting the actual densities: per node and not per face
    std::cout << "Mapping face densities to vertices and applying smoothing..." << std::endl;

    // --- Step 1: Initial area-weighted mapping (same as your original code) ---
    std::vector vertex_densities(n_vertices, 0.0); // Use double for precision
    std::vector ring_areas(n_vertices, 0.0);

    for (size_t i = 0; i < n_faces; ++i) {
        auto& f = faces[i];
        double Ai = face_area(vertices[f.v1], vertices[f.v2], vertices[f.v3]);
        
        vertex_densities[f.v1] += q_arr[i] * Ai;
        ring_areas[f.v1] += Ai;

        vertex_densities[f.v2] += q_arr[i] * Ai;
        ring_areas[f.v2] += Ai;

        vertex_densities[f.v3] += q_arr[i] * Ai;
        ring_areas[f.v3] += Ai;
    }

    for (size_t i = 0; i < n_vertices; ++i) {
        if (ring_areas[i] > 1e-12) { // Avoid division by zero
            vertex_densities[i] /= ring_areas[i];
        }
    }
    ring_areas.clear();
    free(q_arr);



    std::vector<std::vector<int>> adjacency(n_vertices);
    for (const auto& face : faces) {

        // Add edges to adjacency list, avoiding duplicates
        adjacency[face.v1].push_back(face.v2); adjacency[face.v1].push_back(face.v3);
        adjacency[face.v2].push_back(face.v1); adjacency[face.v2].push_back(face.v3);
        adjacency[face.v3].push_back(face.v1); adjacency[face.v3].push_back(face.v2);
    }
    // Clean up duplicates
    for(auto& neighbors : adjacency) {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }



    const int smoothing_iterations = 1; 
    std::vector<double> smoothed_densities = vertex_densities; // Work on a copy

    for (int iter = 0; iter < smoothing_iterations; ++iter) {
        for (size_t i = 0; i < n_vertices; ++i) {
            if (adjacency[i].empty()) continue;

            double neighbor_sum = 0.0;
            for (int neighbor_id : adjacency[i]) {
                neighbor_sum += vertex_densities[neighbor_id];
            }
            smoothed_densities[i] = neighbor_sum / adjacency[i].size();
        }
        vertex_densities = smoothed_densities; // Update for the next iteration
    }


    // --- Final Step: Assign smoothed densities back to the main vertex data structure ---
    for (auto& v : vertices) {
        v.density = vertex_densities[v.id];
    }

    std::cout << "Smoothing complete." << std::endl;



    write_vtu(mesh_file, vertices, faces, n_vertices, n_faces);


*/
