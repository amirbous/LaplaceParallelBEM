# Root CMakeLists.txt

cmake_minimum_required(VERSION 3.10)
project(Par3DLaplaceBEM)

# Set global C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(GINKGO_USR "/storage/home/abouslama/ginkgoParallel/ginkgo")

# --- User-configurable path for Ginkgo
set(Ginkgo_DIR "${GINKGO_USR}/build/install/lib64/cmake/Ginkgo" CACHE PATH "Path to the Ginkgo package config directory")

# --- Add the subdirectories
# The copy command has been moved from here into the subdirectories.
add_subdirectory(serial)
add_subdirectory(parallel)

# --- Add tests
include(CTest) 
add_subdirectory(test)
