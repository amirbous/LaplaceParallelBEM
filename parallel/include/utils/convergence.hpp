/*******************************************************************************
 * 
 * 
 * 
 *  File with functionality to check for convergence
 *  copy paste after main routine in main in case to check for convergence
 *  Cylinder models can be used for this, and ./test_convergence script can 
 *  run directly
 * 
 *  !!!!!! This file is not built or compiled, it is solely a backup for convergence check
 *  whenever that is needed
 * 
 * ********************************************************************************/


    /***************************************************
     *
     *
     *      This part is only to ensure the convergence of solver
     *      will be deleted or left for backup once convergence is correct
     *
     *
     *      
    *****************************************************/
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

