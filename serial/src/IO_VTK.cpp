#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>

#include <iomanip>

#include "geometry.hpp"


void read_vtk(const std::string problem_name, std::vector<Vertex>& vertices, std::vector<Face>& faces,
				int& nvertices, int& nfaces) {


	// used to extract the used faces

	struct Face_vid {
		int v1_id;
		int v2_id;
		int v3_id;
	};

	size_t nall_nodes{0}, ncellsOffsets{0}, 
				global_vertexId{0};

	float xi{0.0f}, yi{0.0f}, zi{0.0f};
	int v1_i{0}, v2_i{0}, v3_i{0}, face_i{0}, face_j{0};
	int offseti{0};


	int val_buffer{0};

	std::ifstream fstream;
	std::string line_buffer, word_buffer;

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


	std::vector<Vertex> all_vertices(nall_nodes);

	// to extract only the main vertices
	char *is_main_vertex = (char *) calloc(nall_nodes, sizeof(char));


	for (size_t v_id = 0; v_id < nall_nodes; v_id++) {

		fstream >> xi;
		fstream >> yi;
		fstream >> zi;
		all_vertices[v_id] = Vertex(xi, yi, zi, v_id);


	}

	// parse and store cells offsets
	fstream >> word_buffer;

	fstream >> nfaces;
	fstream >> ncellsOffsets;

	nfaces = nfaces - 1;

	std::vector<Face_vid> faces_vids(nfaces);


	faces = std::vector<Face>(nfaces);


	fstream >> word_buffer;
	fstream >> word_buffer;


	for (size_t offset_id = 0; offset_id < nfaces + 1; offset_id++) {
		fstream >> offseti;

	}

	fstream >> word_buffer;
	fstream >> word_buffer;


	for (size_t face_id = 0; face_id < nfaces; face_id++ ) {


			fstream >> v1_i;
			fstream >> v2_i;
			fstream >> v3_i;

    		faces_vids[face_id] = Face_vid{v1_i, v2_i, v3_i};

			fstream >> val_buffer;
			fstream >> val_buffer;
			fstream >> val_buffer;

			is_main_vertex[v1_i] = 1;
			is_main_vertex[v2_i] = 1;
			is_main_vertex[v3_i] = 1;
	}



	for (int i = 0; i < nall_nodes; i++) {
		nvertices += is_main_vertex[i];
	}


	vertices = std::vector<Vertex>(nvertices);

	std::map<int, int> old_new_vertexId;

	for (int i = 0; i < nall_nodes; i++) {
    	if (is_main_vertex[i] == 1) {
        	vertices[global_vertexId] = Vertex(all_vertices[i]);
        	old_new_vertexId[i] = global_vertexId++;
    	}
	}
	for (int i = 0; i < nfaces; i++) {
			faces[i] = Face( old_new_vertexId[faces_vids[i].v1_id], 
                   				old_new_vertexId[faces_vids[i].v2_id], 
                   				old_new_vertexId[faces_vids[i].v3_id]
                   			   );
	}
	for (int i = 0; i < nvertices; i++) {
		vertices[i].id = i;
	}


	all_vertices.clear();
	faces_vids.clear();

}


void write_vtu(const std::string problem_name, const std::vector<Vertex>& vertices, const std::vector<Face>& faces,
							int &nvertices, int &nfaces) { 


	std::ofstream fstream;


	std::string fname = problem_name + ".vtu";

	fstream.open(fname);
	if (fstream.is_open()) {

		// headers
		fstream << "<?xml version=\"1.0\"?>" << std::endl;
		fstream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
		fstream << " byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
		fstream << "<UnstructuredGrid>" << std::endl;
		fstream << "<Piece NumberOfPoints=\"" << nvertices << "\"  ";
		fstream << "NumberOfCells=\"" << nfaces << "\">" << std::endl;
		fstream << "<Points>" << std::endl;
		fstream << "<DataArray type=\"Float32\" Name=\"Points\"";
		fstream << " NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;

		// Writing the vertices
		for (size_t i = 0; i < nvertices; i++) {

			fstream << std::fixed << std::setprecision(7) << vertices[i].x 
			 << " " << std::fixed << std::setprecision(7) << vertices[i].y
			 << " " << std::fixed << std::setprecision(7) << vertices[i].z
			 << std::endl;

		}

		fstream << "</DataArray>" << std::endl;
		fstream << "</Points>" << std::endl;
		fstream << "<Cells>" << std::endl;
		fstream << "<DataArray type=\"Int64\" ";
		fstream << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

		for (size_t i = 0; i < nfaces; i++) {
			fstream << faces[i].v1
			 << " " << faces[i].v2
			 << " " << faces[i].v3 
			 << std::endl;
		}

		fstream << "</DataArray>" << std::endl;

		fstream << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
		for (int i = 3; i < nfaces * 3 + 1; i+= 3) {
			fstream << i << ((i % 30) != 0 ? " " : " \n");
		}

		if (nfaces  *  3 % 30 != 0) {
			fstream << std::endl;
		}

		fstream << "</DataArray>" << std::endl;


		fstream << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;

		for (size_t i = 0; i < nfaces; i++) {
			fstream << "5" << ((i+1) % 10 != 0 ? " " : " \n");
		}

		if (nfaces % 10 != 0) {
			fstream << std::endl;
		}
		fstream << "</DataArray>" << std::endl;
		fstream << "</Cells>" << std::endl;

		fstream << "<PointData>" << std::endl;
		fstream << "<DataArray type=\"Float32\" Name=\"Potential\" ";
		fstream << "NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;

		for (int i = 0; i < nvertices; i++) {
			fstream << vertices[i].potential << std::endl;
		}

		fstream << "</DataArray>" << std::endl;

		fstream << "<DataArray type=\"Float32\" Name=\"density\" ";
		fstream << "NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
		for (int i = 0; i < nvertices; i++) {
			fstream << vertices[i].density << std::endl;
		}

		fstream << "</DataArray>" << std::endl;

		fstream << "</PointData>" << std::endl;

		fstream << "</Piece>" << std::endl;
		fstream << "</UnstructuredGrid>" << std::endl;
		fstream << "</VTKFile>";


		fstream.close();


	} else {

		std::cerr << "Error in opening result file" << std::endl;
	}

}

void write_vtu(const std::string problem_name, const std::vector<Vertex>& vertices, const std::vector<Face>& faces,
							int &nvertices, int &nfaces) ;


void write_vtu(const std::string problem_name, const std::vector<Vertex>& vertices, const std::vector<Face>& faces,
							int &nvertices, int &nfaces);


