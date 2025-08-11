cmake_minimum_required(VERSION 3.10)

include(FetchContent)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# --- Copy the shared geometry directory into this target's build folder
file(COPY ${PROJECT_SOURCE_DIR}/geo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

enable_testing()

set(TARGET_GEOMETRY "geometry_test")

add_executable(${TARGET_GEOMETRY}
    ${PROJECT_SOURCE_DIR}/serial/src/IO_VTK.cpp
    ${PROJECT_SOURCE_DIR}/serial/src/geo_math.cpp
    ${PROJECT_SOURCE_DIR}/test/geometry_test.cpp   # Add your test file here
)

target_include_directories(${TARGET_GEOMETRY}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${TARGET_GEOMETRY} gtest_main)

add_test(NAME ${TARGET_GEOMETRY} COMMAND ${TARGET_GEOMETRY})
