#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "../include/geometry.hpp"

template
<typename T>
void read_vtk(const std::string problem_name, Vertex<T>** vertices, Face<T>** faces,
				size_t *nvertices, size_t *nfaces) {

	size_t nall_nodes{0}, ncellsOffsets{0}, 
				global_vertexId{0};



	float xi{0.0f}, yi{0.0f}, zi{0.0f};
	int v1_i{0}, v2_i{0}, v3_i{0}, face_i{0}, face_j{0};
	int offseti{0};


	std::ifstream fstream;
	std::string line_buffer, word_buffer;
	int val_buffer{0};
	std::string fname = "geo/" + problem_name + ".vtk";


	const uint8_t nHeaderLines = 4;


	fstream.open(fname);

	if (!fstream) {
		std::cerr << "Error: file '" << fname << "' could not be opened or does not exist." << std::endl;
		
	}
	
	// skip 3 first lines: useless info
	for (int i = 0; i < nHeaderLines; i++) {
		std::getline(fstream, line_buffer);

	}

	// skip first word
	fstream >> word_buffer;
	fstream >> nall_nodes;
	fstream >> word_buffer;



	Vertex<T>** all_vertices = (Vertex<T>**) malloc(sizeof(Vertex<T>*) * nall_nodes);
	// to extract only the main vertices
	char *is_main_vertex = (char *) calloc(nall_nodes, sizeof(char));


	


	//*vertices = (Vertex<T>*) malloc(sizeof(Vertex<T>) * nvertices);



	for (size_t v_id = 0; v_id < nall_nodes; v_id++) {

		fstream >> xi;
		fstream >> yi;
		fstream >> zi;
		all_vertices[v_id] = new Vertex<T>(v_id, xi, yi, zi);


	}

	// parse and store cells offsets
	fstream >> word_buffer;

	fstream >> *nfaces;
	fstream >> ncellsOffsets;

	*faces = (Face<T>*) malloc(sizeof(Face<T>) * *nfaces);



	fstream >> word_buffer;
	fstream >> word_buffer;


	for (size_t offset_id = 0; offset_id < nfaces; offset_id++) {
		fstream >> offseti;

	}

	fstream >> word_buffer;
	fstream >> word_buffer;


	for (size_t face_id = 0; face_id < nfaces; face_id++ ) {


			fstream >> v1_i;
			fstream >> v2_i;
			fstream >> v3_i;

    		(*faces)[face_id] = Face<T>(all_vertices[v1_i],
                                all_vertices[v2_i],
                                all_vertices[v3_i]);

			fstream >> val_buffer;
			fstream >> val_buffer;
			fstream >> val_buffer;

			is_main_vertex[v1_i] = 1;
			is_main_vertex[v2_i] = 1;
			is_main_vertex[v3_i] = 1;
	}



	for (int i = 0; i < nall_nodes; i++) {
		*nvertices += is_main_vertex[i];
	}


	std::cout << nall_nodes << " - " << nvertices << std::endl;

	*vertices = (Vertex<T>*) malloc(sizeof(Vertex<T>) * *nvertices); 

	for (int i = 0; i < nall_nodes; i++) {
    	if (is_main_vertex[i] == 1) {
        	(*vertices)[global_vertexId] = *all_vertices[i]; 
        	(*vertices)[global_vertexId].id = global_vertexId++;
    	}
	}

	// freeing temporary arrays
	//free(all_vertices);
	//free(is_main_vertex);
}



template void read_vtk<double>(const std::string problem_name, Vertex<double>**, Face<double>**);

