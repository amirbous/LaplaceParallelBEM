cmake_minimum_required(VERSION 3.10)


project(Par3DLaplaceBEM)


set(CMAKE_CXX_STANDARD 17)



# Ginkgo installation: important user specifies
set(GINKGO_DIR "$ENV{HOME}/ginkgo/build/install")

#serial target
set(SERIAL_TARGET "electromain")


include_directories(include)


file(COPY "${PROJECT_SOURCE_DIR}/geo" DESTINATION "${CMAKE_BINARY_DIR}")

find_package(OpenMP REQUIRED)
find_package(Ginkgo REQUIRED PATHS "${GINKGO_DIR}/lib/cmake/ginkgo")



# Sroucese for default version 
SET(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
SET(SOURCES ${SOURCE_DIR}/main.cpp ${SOURCE_DIR}/IO_VTK.cpp ${SOURCE_DIR}/geo_math.cpp)


add_executable(${SERIAL_TARGET} ${SOURCES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(${SERIAL_TARGET} PUBLIC OpenMP::OpenMP_CXX)
endif()

if (Ginkgo_FOUND)
    target_link_libraries(${SERIAL_TARGET} PUBLIC Ginkgo::ginkgo)
else()
    message(WARNING "Ginkgo not found for 'main' target. Proceeding without it.")
endif()


####################################################################
# Setup for MPI
####################################################################


find_package(MPI REQUIRED)


SET(PARALLEL_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src_mpi")
SET(PARALLEL_SOURCES ${PARALLEL_SOURCE_DIR}/main.cpp ${PARALLEL_SOURCE_DIR}/IO_VTK.cpp ${PARALLEL_SOURCE_DIR}/geo_math.cpp)


#serial target
set(PARAL_TARGET "electroparal")

add_executable(${PARAL_TARGET} ${PARALLEL_SOURCES})

# Link libraries for the parallel version, including MPI
target_link_libraries(${PARAL_TARGET} PUBLIC
    MPI::MPI_CXX
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${PARAL_TARGET} PUBLIC OpenMP::OpenMP_CXX)
endif()

if (Ginkgo_FOUND)
    target_link_libraries(${PARAL_TARGET} PUBLIC Ginkgo::ginkgo)
else()
    message(WARNING "Ginkgo not found for 'paral' target. Proceeding without it.")
endif()
