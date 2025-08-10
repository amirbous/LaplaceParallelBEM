cmake_minimum_required(VERSION 3.10)

include(FetchContent)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(geometry_test geometry_test.cpp)
target_link_libraries(geometry_test gtest_main)

add_test(NAME geometry_test COMMAND geometry_test)
