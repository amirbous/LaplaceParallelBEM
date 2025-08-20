/***********************************************
* File for verifying corret data and array correctness
************************************************/

#include <gtest/gtest.h>

#include "IO_VTK.hpp"
#include "geometry.hpp"

#include <iostream>
#include <vector>
#include <string>

TEST(ArrayAssemblyTest, parsing_elementsNumber) {

    int SPHERE_N_FACES = 340;
    int SPHERE_N_VERTICES = 172;
    std::string SPHERE_PROBLEMNAME = "sphere";


    int n_vertices{0}, n_faces{0};
    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    read_vtk(SPHERE_PROBLEMNAME, vertices, faces, n_vertices, n_faces);

    // TODO: faces array and vertices array read from a file
}

TEST(ArrayAssemblyTest, matrixArray_serial) {

    int SPHERE_N_FACES = 340;
    int SPHERE_N_VERTICES = 172;
    std::string SPHERE_PROBLEMNAME = "sphere";

    int n_vertices{0}, n_faces{0};
    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    read_vtk(SPHERE_PROBLEMNAME, vertices, faces, n_vertices, n_faces);

    std::string groundthruth_matrix = "sphere_correct.mtx";

    // TODO: complete matrix flow, then compare to file
    

}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
