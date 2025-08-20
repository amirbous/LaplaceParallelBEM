#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>


#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <cctype>

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

void read_vtu_sol(const std::string problem_name, std::vector<Vertex>& vertices,
                  std::vector<Face>& faces,
                  int& nvertices, int& nfaces) {

    struct Face_vid {
        int v1_id;
        int v2_id;
        int v3_id;
    };

    std::string fname = "geo/" + problem_name + ".vtu";
    std::ifstream fstream(fname);
    if (!fstream) {
        std::cerr << "Error: file '" << fname << "' could not be opened or does not exist." << std::endl;
        return;
    }

    // Read whole file into string for easy substring extraction
    std::ostringstream ss;
    ss << fstream.rdbuf();
    std::string file = ss.str();

    // Helper lambda: find the inner text between an opening tag and its corresponding closing tag.
    auto extract_between = [&](const std::string& open_tag, const std::string& close_tag, size_t start_pos = 0) -> std::string {
        size_t a = file.find(open_tag, start_pos);
        if (a == std::string::npos) return "";
        a = file.find('>', a);
        if (a == std::string::npos) return "";
        a++; // position after '>'
        size_t b = file.find(close_tag, a);
        if (b == std::string::npos) return "";
        return file.substr(a, b - a);
    };

    // --- 1) Points ---
    // Find the <Points>...</Points> block and within it the first DataArray (assumed to be Points)
    size_t points_block_start = file.find("<Points");
    if (points_block_start == std::string::npos) {
        std::cerr << "Error: <Points> block not found in '" << fname << "'." << std::endl;
        return;
    }
    std::string points_data = extract_between("<Points", "</Points>", points_block_start);
    if (points_data.empty()) {
        std::cerr << "Error: could not extract Points DataArray." << std::endl;
        return;
    }

    // Within points_data, find the DataArray section (there may be attributes on the opening tag)
    size_t da_start = points_data.find("<DataArray");
    if (da_start == std::string::npos) {
        std::cerr << "Error: Points DataArray not found." << std::endl;
        return;
    }
    // We want content between that DataArray's '>' and its closing tag
    size_t da_real_start = points_data.find('>', da_start);
    size_t da_end = points_data.find("</DataArray>", da_real_start);
    if (da_real_start == std::string::npos || da_end == std::string::npos) {
        std::cerr << "Error: Points DataArray malformed." << std::endl;
        return;
    }
    std::string points_inner = points_data.substr(da_real_start + 1, da_end - (da_real_start + 1));

    // Parse floats triplets
    std::istringstream pts_stream(points_inner);
    float xi = 0.0f, yi = 0.0f, zi = 0.0f;
    std::vector<Vertex> all_vertices;
    int temp_id = 0;
    while (pts_stream >> xi >> yi >> zi) {
        all_vertices.emplace_back(xi, yi, zi, temp_id++);
    }
    size_t nall_nodes = all_vertices.size();
    if (nall_nodes == 0) {
        std::cerr << "Warning: no points read from '" << fname << "'." << std::endl;
    }

    // --- 2) Cells: connectivity ---
    size_t cells_block_start = file.find("<Cells");
    if (cells_block_start == std::string::npos) {
        std::cerr << "Error: <Cells> block not found in '" << fname << "'." << std::endl;
        return;
    }
    std::string cells_block = extract_between("<Cells", "</Cells>", cells_block_start);
    if (cells_block.empty()) {
        std::cerr << "Error: could not extract Cells block." << std::endl;
        return;
    }

    // Find connectivity DataArray in cells block
    size_t conn_pos = cells_block.find("Name=\"connectivity\"");
    if (conn_pos == std::string::npos) {
        // try without Name attr (just take first DataArray inside Cells)
        conn_pos = cells_block.find("<DataArray");
        if (conn_pos == std::string::npos) {
            std::cerr << "Error: connectivity DataArray not found in Cells." << std::endl;
            return;
        }
    } else {
        // move backwards to the opening <DataArray
        conn_pos = cells_block.rfind("<DataArray", conn_pos);
    }

    size_t conn_da_start = cells_block.find('>', conn_pos);
    size_t conn_da_end = cells_block.find("</DataArray>", conn_da_start);
    if (conn_da_start == std::string::npos || conn_da_end == std::string::npos) {
        std::cerr << "Error: Cells connectivity DataArray malformed." << std::endl;
        return;
    }
    std::string conn_inner = cells_block.substr(conn_da_start + 1, conn_da_end - (conn_da_start + 1));

    // Parse integers triples (connectivity)
    std::istringstream conn_stream(conn_inner);
    std::vector<Face_vid> faces_vids;
    int a = 0, b = 0, c = 0;
    while (conn_stream >> a >> b >> c) {
        faces_vids.push_back(Face_vid{a, b, c});
    }

    nfaces = static_cast<int>(faces_vids.size());
    if (nfaces == 0) {
        std::cerr << "Warning: no faces read from '" << fname << "'." << std::endl;
    }

    // --- 3) PointData: find Electric Field Array ---
    // Find the <PointData ...> block (if present), then search for DataArray whose Name contains "Electric" (case-insensitive)
    size_t pointdata_pos = file.find("<PointData");
    std::vector<float> electric_field_values;
    if (pointdata_pos != std::string::npos) {
        // search for DataArray occurrences within PointData block
        size_t pd_end = file.find("</PointData>", pointdata_pos);
        if (pd_end == std::string::npos) pd_end = pointdata_pos + 2000; // fallback
        std::string pd_block = file.substr(pointdata_pos, pd_end - pointdata_pos);

        // find DataArray that contains "Electric" in its opening tag
        size_t search_pos = 0;
        size_t found_da = std::string::npos;
        while (true) {
            size_t da = pd_block.find("<DataArray", search_pos);
            if (da == std::string::npos) break;
            size_t da_close = pd_block.find('>', da);
            if (da_close == std::string::npos) break;
            std::string opening = pd_block.substr(da, da_close - da + 1);
            // lowercase copy for case-insensitive search
            std::string opening_l = opening;
            std::transform(opening_l.begin(), opening_l.end(), opening_l.begin(),
                           [](unsigned char ch) { return std::tolower(ch); });
            if (opening_l.find("electric") != std::string::npos) {
                found_da = da;
                break;
            }
            search_pos = da_close + 1;
        }

        if (found_da != std::string::npos) {
            size_t ef_da_start = pd_block.find('>', found_da);
            size_t ef_da_end = pd_block.find("</DataArray>", ef_da_start);
            if (ef_da_start != std::string::npos && ef_da_end != std::string::npos) {
                std::string ef_inner = pd_block.substr(ef_da_start + 1, ef_da_end - (ef_da_start + 1));
                std::istringstream ef_stream(ef_inner);
                float val = 0.0f;
                while (ef_stream >> val) electric_field_values.push_back(val);
            }
        } else {
            // fallback: maybe named '|Electric Field|' with extra spaces, try a substring search in pd_block
            size_t alt = pd_block.find("Electric");
            if (alt != std::string::npos) {
                // extract the whole DataArray surrounding it
                size_t da_open = pd_block.rfind("<DataArray", alt);
                size_t da_close = pd_block.find("</DataArray>", alt);
                if (da_open != std::string::npos && da_close != std::string::npos) {
                    size_t da_real_start = pd_block.find('>', da_open);
                    if (da_real_start != std::string::npos) {
                        std::string ef_inner = pd_block.substr(da_real_start + 1, da_close - (da_real_start + 1));
                        std::istringstream ef_stream(ef_inner);
                        float val = 0.0f;
                        while (ef_stream >> val) electric_field_values.push_back(val);
                    }
                }
            }
        }
    } // else: no PointData found

    // If we failed to read electric field values, fill with zeros to avoid out-of-range.
    if (electric_field_values.size() < nall_nodes) {
        electric_field_values.resize(nall_nodes, 0.0f);
    }

    // --- 4) Determine main vertices (those used by faces) ---
    std::vector<char> is_main_vertex(nall_nodes, 0);
    for (const auto &fv : faces_vids) {
        if (fv.v1_id >= 0 && fv.v1_id < (int)nall_nodes) is_main_vertex[fv.v1_id] = 1;
        if (fv.v2_id >= 0 && fv.v2_id < (int)nall_nodes) is_main_vertex[fv.v2_id] = 1;
        if (fv.v3_id >= 0 && fv.v3_id < (int)nall_nodes) is_main_vertex[fv.v3_id] = 1;
    }

    // Count vertices that are flagged
    nvertices = 0;
    for (size_t i = 0; i < nall_nodes; ++i) nvertices += is_main_vertex[i];

    // Build new vertices array and map old -> new indices
    vertices = std::vector<Vertex>(nvertices);
    std::map<int,int> old_new_vertexId;
    int global_vertexId = 0;
    for (int i = 0; i < (int)nall_nodes; ++i) {
        if (is_main_vertex[i]) {
            vertices[global_vertexId] = Vertex(all_vertices[i]); // copy coordinates
            // assign electric field value to vertex attribute
            // NOTE: adjust member name if your Vertex uses a different name than 'densities'
            vertices[global_vertexId].density = electric_field_values[i];
            old_new_vertexId[i] = global_vertexId++;
        }
    }

    // Build faces with remapped vertex ids
    faces = std::vector<Face>(nfaces);
    for (int i = 0; i < nfaces; ++i) {
        int a_old = faces_vids[i].v1_id;
        int b_old = faces_vids[i].v2_id;
        int c_old = faces_vids[i].v3_id;
        int a_new = old_new_vertexId.count(a_old) ? old_new_vertexId[a_old] : -1;
        int b_new = old_new_vertexId.count(b_old) ? old_new_vertexId[b_old] : -1;
        int c_new = old_new_vertexId.count(c_old) ? old_new_vertexId[c_old] : -1;
        if (a_new < 0 || b_new < 0 || c_new < 0) {
            // skip or create an invalid face; here we create with zeros to preserve indexing
            faces[i] = Face(0,0,0);
        } else {
            faces[i] = Face(a_new, b_new, c_new);
        }
    }

    // Reset vertex ids to contiguous 0..nvertices-1
    for (int i = 0; i < nvertices; ++i) {
        vertices[i].id = i;
    }

    // done
    // Note: all_vertices and faces_vids will be destroyed on function exit
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


