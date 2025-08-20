/********************************************************
* File for testing correctness of mathematical operators
********************************************************/
#include <gtest/gtest.h>

#include "IO_VTK.hpp"
#include "geometry.hpp"
#include "geo_math.hpp"

#include <iostream>
#include <vector>
#include <string>

TEST(NumericTest, area_computation) {
    Vertex v1 = Vertex(1, 2, 3);
    Vertex v2 = Vertex(-3, 4, 5);
    Vertex v3 = Vertex(6, -3, 4);


    float compute_area = face_area(v1, v2, v3);

    // computed by hand then rounded in python
    float const true_area{10.488089};

    ASSERT_FLOAT_EQ(compute_area, true_area) << "computed area " << compute_area << " different from true area " << true_area;

}
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
