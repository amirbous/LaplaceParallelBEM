cmake_minimum_required(VERSION 3.10)

# project name
project(Par3DLaplaceBEM)

# c++ standard
set(CMAKE_CXX_STANDARD 17)

set(GINKGO_DIR "$ENV{HOME}/ginkgo/build/install")


# sources to compile
SET(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")

# sources to compile
SET(SOURCES ${SOURCE_DIR}/main.cpp ${SOURCE_DIR}/IO_VTK.cpp ${SOURCE_DIR}/geo_math.cpp)

# creates a directory for geometries
file(COPY "${PROJECT_SOURCE_DIR}/geo" DESTINATION "${CMAKE_BINARY_DIR}")

add_executable(main ${SOURCES})

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)

find_package(Ginkgo REQUIRED PATHS "${GINKGO_DIR}/lib/cmake/ginkgo")


if(OpenMP_CXX_FOUND)
    target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)
endif()

if (Ginkgo_FOUND)
    target_link_libraries(main PUBLIC Ginkgo::ginkgo)
else()
    message(WARNING "Ginkgo not found. Proceeding without it.")
endif()

target_link_libraries(main PUBLIC Eigen3::Eigen)

# include directoy for hpp files
include_directories(include)
