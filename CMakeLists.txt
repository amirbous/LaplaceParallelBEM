# Root CMakeLists.txt

cmake_minimum_required(VERSION 3.10)
project(Par3DLaplaceBEM)

# Set global C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- User-configurable path for Ginkgo
set(GINKGO_DIR "$ENV{HOME}/ginkgo/build/install" CACHE PATH "Path to the Ginkgo installation directory")

# --- Add the subdirectories
# The copy command has been moved from here into the subdirectories.
add_subdirectory(serial)
add_subdirectory(parallel)

# --- Add tests
include(CTest) 
add_subdirectory(test)