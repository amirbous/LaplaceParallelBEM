#include <iostream>
#include <vector>
#include <string>

#include <chrono>
#include <limits>


#include <ginkgo/ginkgo.hpp>

#include "../include/geometry.hpp"
#include "../include/IO_VTK.hpp"
#include "../include/geo_math.hpp"


int main(int argc, char* argv[]) {
    using ValueType = double;
    std::string mesh_file = (argc > 1 ? argv[1] : "sphere");

    size_t n_vertices{0}, n_faces{0};
    std::vector<Vertex<ValueType>> vertices;
    std::vector<Face<ValueType>> faces;



    float mesh_centroid[] = {0.0f, 0.0f, 0.0f};
    float cent[] = {0.0f, 0.0f, 0.0f};



    read_vtk<ValueType>(mesh_file, vertices, faces, n_vertices, n_faces);
    std::cout << "Mesh: " << mesh_file << ", " << n_vertices << " nodes, " << n_faces << " faces" << std::endl;



    ValueType* G_arr = (ValueType *) calloc(n_faces * n_faces, sizeof(ValueType));
    ValueType* phi_arr = (ValueType *) calloc(n_faces, sizeof(ValueType));
    ValueType* q_arr = (ValueType *) calloc(n_faces, sizeof(ValueType));


    // Boundary Conditions
    //for (auto& v : vertices) { v.potential = 10;} 
    for (auto& v : vertices) { v.potential = 10;} 
    
    // compute the center of the object


    std::cout << "Started assembling the matrix" << std::endl;
    
    /*****************************************************************************
    *
    *
    *             main loop to assemble the matrix: TODO parallelize
    *
    *
    ******************************************************************************/


    std::chrono::steady_clock::time_point begin_assembleMatrix = std::chrono::steady_clock::now(); // measure time for this section


    for (size_t i = 0; i < n_faces; ++i) {
        cent[0] = (faces[i].v1->x + faces[i].v2->x + faces[i].v3->x) / 3;
        cent[1] = (faces[i].v1->y + faces[i].v2->y + faces[i].v3->y) / 3;
        cent[2] = (faces[i].v1->z + faces[i].v2->z + faces[i].v3->z) / 3;

        for (size_t j = 0; j < n_faces; ++j) {
            if (i == j) {
                // compute area of face:
                G_arr[n_faces * i + i] = regularized_integral(faces[j].v1, faces[j].v2, faces[j].v3);

            } else {
                // evaluating the gaussian integral
                G_arr[n_faces * i + j] = gauss_integral(faces[j].v1, faces[j].v2, faces[j].v3, cent);
            }
        }
    }

    std::chrono::steady_clock::time_point end_assembleMatrix = std::chrono::steady_clock::now(); // end

    /******************************************************************************
    ******************************************************************************/


    auto time_toAssemble = std::chrono::duration_cast<std::chrono::milliseconds>(end_assembleMatrix - begin_assembleMatrix).count();
    std::cout << "finished assembling the matrix in: " << time_toAssemble << "ms" << std::endl;

    // constructing load vector ==> density averaged over faces
    for (size_t i = 0; i < n_faces; ++i) {
        phi_arr[i] = (faces[i].v1->potential + faces[i].v2->potential + faces[i].v3->potential) / 3.0;
    }


    // solving the system
    /*--------------------------------------------------
            Solving using ginkgo
    ---------------------------------------------------*/


    std::cout << gko::version_info::get() << std::endl;
    const auto exec = gko::OmpExecutor::create();

    //system matrix
    std::shared_ptr<gko::matrix::Dense<ValueType>> G = 
            gko::matrix::Dense<ValueType>::create(                             
                exec, gko::dim<2>{n_faces, n_faces},      
                gko::array<ValueType>::view(exec, n_faces * n_faces, G_arr)
                , n_faces); 

    //RHS, load vector
    std::shared_ptr<gko::matrix::Dense<ValueType>> phi = 
            gko::matrix::Dense<ValueType>::create(                             
                exec, gko::dim<2>{n_faces, 1},      
                gko::array<ValueType>::view(exec, n_faces, phi_arr)
                , 1); 

    // to store solution
    std::shared_ptr<gko::matrix::Dense<ValueType>> q = 
            gko::matrix::Dense<ValueType>::create(                             
                exec, gko::dim<2>{n_faces, 1},      
                gko::array<ValueType>::view(exec, n_faces, q_arr)
                , 1);

    auto gmres_gen = gko::solver::Gmres<ValueType>::build()
        .with_criteria(
            gko::stop::Iteration::build().with_max_iters(200u).on(exec),
            gko::stop::ResidualNorm<ValueType>::build()
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
    std::vector<double> vertex_densities(n_vertices, 0.0); // Use double for precision
    std::vector<double> ring_areas(n_vertices, 0.0);

    for (size_t i = 0; i < n_faces; ++i) {
        auto& f = faces[i];
        double Ai = face_area(f.v1, f.v2, f.v3);
        
        vertex_densities[f.v1->id] += q_arr[i] * Ai;
        ring_areas[f.v1->id] += Ai;

        vertex_densities[f.v2->id] += q_arr[i] * Ai;
        ring_areas[f.v2->id] += Ai;

        vertex_densities[f.v3->id] += q_arr[i] * Ai;
        ring_areas[f.v3->id] += Ai;
    }

    for (size_t i = 0; i < n_vertices; ++i) {
        if (ring_areas[i] > 1e-12) { // Avoid division by zero
            vertex_densities[i] /= ring_areas[i];
        }
    }
    ring_areas.clear();
    free(q_arr);


    // --- Step 2: Build vertex adjacency list for neighbor finding ---
    std::vector<std::vector<int>> adjacency(n_vertices);
    for (const auto& face : faces) {
        int v1_id = face.v1->id;
        int v2_id = face.v2->id;
        int v3_id = face.v3->id;
        // Add edges to adjacency list, avoiding duplicates
        adjacency[v1_id].push_back(v2_id); adjacency[v1_id].push_back(v3_id);
        adjacency[v2_id].push_back(v1_id); adjacency[v2_id].push_back(v3_id);
        adjacency[v3_id].push_back(v1_id); adjacency[v3_id].push_back(v2_id);
    }
    // Clean up duplicates
    for(auto& neighbors : adjacency) {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }


    // --- Step 3: Apply Laplacian Smoothing ---
    const int smoothing_iterations = 2; // Tune this parameter (5-20 is a good range)
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



    write_vtu<ValueType>(mesh_file, vertices, faces, n_vertices, n_faces);

    size_t n_vertices_verif{0}, n_faces_verif{0};
    std::vector<Vertex<ValueType>> vertices_verif;
    std::vector<Face<ValueType>> faces_verif;

    std::string mesh_file_verif = "Cylinder_7";


    read_vtk<ValueType>(mesh_file_verif, vertices_verif, faces_verif, n_vertices_verif, n_faces_verif);

    ValueType* G_arr_verif = (ValueType *) calloc(n_faces_verif * n_faces_verif, sizeof(ValueType));
    ValueType* phi_arr_verif = (ValueType *) calloc(n_faces_verif, sizeof(ValueType));
    ValueType* q_arr_verif = (ValueType *) calloc(n_faces_verif, sizeof(ValueType));


    // Boundary Conditions
    //for (auto& v : vertices) { v.potential = 10;} 
    for (auto& v : vertices_verif) { v.potential = 10;} 
    


    std::cout << "Started assembling verification matrix" << std::endl;
    

    for (size_t i = 0; i < n_faces_verif; ++i) {
        cent[0] = (faces_verif[i].v1->x + faces_verif[i].v2->x + faces_verif[i].v3->x) / 3;
        cent[1] = (faces_verif[i].v1->y + faces_verif[i].v2->y + faces_verif[i].v3->y) / 3;
        cent[2] = (faces_verif[i].v1->z + faces_verif[i].v2->z + faces_verif[i].v3->z) / 3;

        for (size_t j = 0; j < n_faces_verif; ++j) {
            if (i == j) {
                // compute area of face:
                G_arr_verif[n_faces_verif * i + i] = regularized_integral(faces_verif[j].v1, faces_verif[j].v2, faces_verif[j].v3);

            } else {
                // evaluating the gaussian integral
                G_arr_verif[n_faces_verif * i + j] = gauss_integral(faces_verif[j].v1, faces_verif[j].v2, faces_verif[j].v3, cent);
            }
        }
    }



    std::cout << "finished assembling verification matrix" << std::endl;

    // constructing load vector ==> density averaged over faces
    for (size_t i = 0; i < n_faces_verif; ++i) {
        phi_arr_verif[i] = (faces_verif[i].v1->potential + faces_verif[i].v2->potential + faces_verif[i].v3->potential) / 3.0;
    }

    //system matrix
    std::shared_ptr<gko::matrix::Dense<ValueType>> G_verif = 
            gko::matrix::Dense<ValueType>::create(                             
                exec, gko::dim<2>{n_faces_verif, n_faces_verif},      
                gko::array<ValueType>::view(exec, n_faces_verif * n_faces_verif, G_arr_verif)
                , n_faces_verif); 

    //RHS, load vector
    std::shared_ptr<gko::matrix::Dense<ValueType>> phi_verif = 
            gko::matrix::Dense<ValueType>::create(                             
                exec, gko::dim<2>{n_faces_verif, 1},      
                gko::array<ValueType>::view(exec, n_faces_verif, phi_arr_verif)
                , 1); 

    // to store solution
    std::shared_ptr<gko::matrix::Dense<ValueType>> q_verif = 
            gko::matrix::Dense<ValueType>::create(                             
                exec, gko::dim<2>{n_faces_verif, 1},      
                gko::array<ValueType>::view(exec, n_faces_verif, q_arr_verif)
                , 1);

    auto gmres_solver_verif = gmres_gen->generate(G_verif);

    gmres_solver_verif->apply(phi_verif, q_verif);

    std::cout << "finished solving vertification matrix" << std::endl;

    free(G_arr_verif);
    free(phi_arr_verif);

    // --- Step 1: Initial area-weighted mapping (same as your original code) ---
    std::vector<double> vertex_densities_verif(n_vertices_verif, 0.0); // Use double for precision
    std::vector<double> ring_areas_verif(n_vertices_verif, 0.0);

    for (size_t i = 0; i < n_faces_verif; ++i) {
        auto& f = faces_verif[i];
        double Ai = face_area(f.v1, f.v2, f.v3);
        
        vertex_densities_verif[f.v1->id] += q_arr_verif[i] * Ai;
        ring_areas_verif[f.v1->id] += Ai;

        vertex_densities_verif[f.v2->id] += q_arr_verif[i] * Ai;
        ring_areas_verif[f.v2->id] += Ai;

        vertex_densities_verif[f.v3->id] += q_arr_verif[i] * Ai;
        ring_areas_verif[f.v3->id] += Ai;
    }

    for (size_t i = 0; i < n_vertices_verif; ++i) {
        if (ring_areas_verif[i] > 1e-12) { // Avoid division by zero
            vertex_densities_verif[i] /= ring_areas_verif[i];
        }
    }
    ring_areas_verif.clear();
    free(q_arr_verif);


    // --- Step 2: Build vertex adjacency list for neighbor finding ---
    std::vector<std::vector<int>> adjacency_verif(n_vertices_verif);
    for (const auto& face : faces_verif) {
        int v1_id = face.v1->id;
        int v2_id = face.v2->id;
        int v3_id = face.v3->id;
        // Add edges to adjacency list, avoiding duplicates
        adjacency_verif[v1_id].push_back(v2_id); adjacency_verif[v1_id].push_back(v3_id);
        adjacency_verif[v2_id].push_back(v1_id); adjacency_verif[v2_id].push_back(v3_id);
        adjacency_verif[v3_id].push_back(v1_id); adjacency_verif[v3_id].push_back(v2_id);
    }
    // Clean up duplicates
    for(auto& neighbors : adjacency_verif) {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }


    // --- Step 3: Apply Laplacian Smoothing ---
    std::vector<double> smoothed_densities_verif = vertex_densities_verif; // Work on a copy

    for (int iter = 0; iter < smoothing_iterations; ++iter) {
        for (size_t i = 0; i < n_vertices_verif; ++i) {
            if (adjacency_verif[i].empty()) continue;

            double neighbor_sum = 0.0;
            for (int neighbor_id : adjacency_verif[i]) {
                neighbor_sum += vertex_densities_verif[neighbor_id];
            }
            smoothed_densities_verif[i] = neighbor_sum / adjacency_verif[i].size();
        }
        vertex_densities_verif = smoothed_densities_verif; // Update for the next iteration
    }


    // --- Final Step: Assign smoothed densities back to the main vertex data structure ---
    for (auto& v : vertices_verif) {
        v.density = vertex_densities_verif[v.id];
    }

    std::cout << "Smoothing complete." << std::endl;



    // evaluting the error
    double err{0.0};
    double min_distance{0.0}, mesh_distance{0.0}; // for getting the minimum distance
    int close_v_id{0};


    for (size_t v_id = 0; v_id < n_vertices; v_id++) {

        min_distance = std::numeric_limits<double>::infinity();
        close_v_id = 0;
        for (size_t verif_v_id = 0; verif_v_id < n_vertices_verif; verif_v_id++) {
            
            mesh_distance = std::sqrt(std::pow(vertices[v_id].x - vertices_verif[verif_v_id].x, 2) + 
                                        std::pow(vertices[v_id].y - vertices_verif[verif_v_id].y, 2) +
                                        std::pow(vertices[v_id].z - vertices_verif[verif_v_id].z, 2));
            if (mesh_distance <= min_distance) {
                min_distance = mesh_distance;
                close_v_id = verif_v_id;
            }
        }
        err += std::pow(vertices[v_id].density - vertices_verif[close_v_id].density, 2);

    }
    err = std::sqrt(err / n_vertices);
    std::cout << "RMSE = " << err << std::endl;



    return 0;
}