# parallel/CMakeLists.txt

# --- Find required packages for the parallel version
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Ginkgo REQUIRED)

set(MPI_BENCHMARK_SCRIPT "ElectroMPIJob.sh")
set(HYBRID_BENCHMARK_SCRIPT "ElectroHybridJob.sh")

set(GEOMETRY_DIR ${PROJECT_SOURCE_DIR}/geo)

set(TARGET_NAME "electroparal")
# --- Copy the shared geometry directory into this target's build folder
file(COPY ${PROJECT_SOURCE_DIR}/geo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# --- Define the executable and its source files
add_executable(${TARGET_NAME}
    src_mpi/main.cpp
    src_mpi/IO_VTK.cpp
    src_mpi/geo_math.cpp
)

# --- Add include directories for the target
target_include_directories(${TARGET_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}
)


# --- Link libraries to the target
target_link_libraries(${TARGET_NAME} PUBLIC
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
    Ginkgo::ginkgo
)


file(COPY ${MPI_BENCHMARK_SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${HYBRID_BENCHMARK_SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
