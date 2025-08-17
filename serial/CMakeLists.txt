# serial/CMakeLists.txt

# --- Find required packages for the serial version
find_package(OpenMP REQUIRED)
find_package(Ginkgo REQUIRED)


set(BENCHMARK_SCRIPT "ElectroOMPJob.sh")
set(GEOMETRY_DIR ${PROJECT_SOURCE_DIR}/geo)
set(TARGET_NAME "electromain")

# --- Copy the shared geometry directory into this target's build folder
file(COPY ${GEOMETRY_DIR} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# --- Define the executable and its source files
add_executable(${TARGET_NAME}
    src/main.cpp
    src/IO_VTK.cpp
    src/geo_math.cpp
)

# --- Add include directories for the target
target_include_directories(${TARGET_NAME} PRIVATE ../include)

# --- Link libraries to the target
target_link_libraries(${TARGET_NAME} PUBLIC
    OpenMP::OpenMP_CXX
    Ginkgo::ginkgo
)

# --- Copy benchmarking script
file(COPY ${BENCHMARK_SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
