#  This file is a part of HAXX
#
#  Copyright (c) 2017 David Williams-Young
#  All rights reserved.
#
#  See LICENSE.txt 
#
cmake_minimum_required(VERSION 3.0)
project(HAXX CXX C Fortran)
set(HAXX_VERSION_MAJOR ALPHA)
set(HAXX_VERSION_MINOR 1)
set(HAXX_VERSION_PATCH 0)


include_directories("${PROJECT_SOURCE_DIR}/include")
#include_directories("${PROJECT_BINARY_DIR}/include")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  message(STATUS "${CMAKE_CXX_COMPILER_ID} Recognized as CLANG: Adding -D_HAXX_UT_BUTF_NINCLUDED")
  add_definitions("-D_HAXX_UT_BUTF_NINCLUDED")
endif()

include(HAXXCompilerFlags)

option( HAXX_ENABLE_BENCHMARK "Enable build of HAXX / HBLAS Benchmark" OFF )
option( HAXX_ENABLE_TESTING   "Enable build of HAXX / HBLAS Unit Tests" ON )


add_subdirectory(src)

if( HAXX_ENABLE_TESTING )
  enable_testing()
  add_subdirectory(tests)
endif()


# Install Headers
install(FILES ${PROJECT_SOURCE_DIR}/include/haxx.hpp DESTINATION include)
install(FILES ${PROJECT_SOURCE_DIR}/include/hblas.hpp DESTINATION include)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/haxx DESTINATION include)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/hblas DESTINATION include)
