cmake_minimum_required(VERSION 3.10)

# project name
project(Par3DLaplaceBEM)

# c++ standard
set(CMAKE_CXX_STANDARD 17)

# sources to compile
SET(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")

# sources to compile
SET(SOURCES ${SOURCE_DIR}/main.cpp ${SOURCE_DIR}/IO_VTK.cpp)

# creates a directory for geometries
file(COPY "${PROJECT_SOURCE_DIR}/geo" DESTINATION "${CMAKE_BINARY_DIR}")

# linking to one file
# execution and usage guidelines will be specified
add_executable(main ${SOURCES})

# include directoy for hpp files
include_directories(include)
