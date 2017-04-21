#  This file is a part of HAXX
#
#  Copyright (c) 2017 David Williams-Young
#  All rights reserved.
#
#  See LICENSE.txt 
#
cmake_minimum_required(VERSION 3.0)
project(HAXX CXX C)
set(HAXX_VERSION_MAJOR ALPHA)
set(HAXX_VERSION_MINOR 1)
set(HAXX_VERSION_PATCH 0)


include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  message(STATUS "${CMAKE_CXX_COMPILER_ID} Recognized as CLANG: Adding -D_HAXX_UT_BUTF_NINCLUDED")
  add_definitions("-D_HAXX_UT_BUTF_NINCLUDED")
endif()

enable_testing()
add_subdirectory(src)
