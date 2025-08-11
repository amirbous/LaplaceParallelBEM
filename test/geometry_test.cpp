#include <gtest/gtest.h>

#include "IO_VTK.hpp"
#include "geometry.hpp"

#include <iostream>
#include <vector>
#include <string>

TEST(GeometryTest, parsing_elementsNumber) {

    int SPHERE_N_FACES = 340;
    int SPHERE_N_VERTICES = 172;
    std::string SPHERE_PROBLEMNAME = "sphere";


    int n_vertices{0}, n_faces{0};
    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    read_vtk(SPHERE_PROBLEMNAME, vertices, faces, n_vertices, n_faces);

    ASSERT_EQ(n_vertices, SPHERE_N_VERTICES);
    ASSERT_EQ(n_faces, SPHERE_N_FACES);
}



int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
