cmake_minimum_required(VERSION 3.0)
project(HAXX CXX C)
set(HAXX_VERSION_MAJOR ALPHA)
set(HAXX_VERSION_MINOR 1)
set(HAXX_VERSION_PATCH 0)


include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")

enable_testing()
add_subdirectory(src)
