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

    int n_vertices_verif{0}, n_faces_verif{0};
    std::vector<Vertex> vertices_verif;
    std::vector<Face> faces_verif;

    read_vtu_sol("Trax", vertices_verif, faces_verif, n_vertices_verif, n_faces_verif);
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


    write_vtu("Trax", vertices_verif, faces_verif, n_vertices_verif, n_faces_verif);
