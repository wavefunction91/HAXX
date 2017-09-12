#
#  This file is a part of HAXX
#
#  Copyright (c) 2017 David Williams-Young
#  All rights reserved.
#
#  See LICENSE.txt 
#
#set(Boost_USE_STATIC_LIBS ON)


message( STATUS "HAXX WILL BUILD BENCHMARK EXECUTABLES" )

add_executable(gemm gemm.cxx)

set( DOT_SRC dot.cxx )

# Make sure that the fortran functions get compiled for
# the benchmark
if( NOT ENABLE_GENERIC_FORTRAN )
  set( DOT_SRC ${DOT_SRC} 
    ${PROJECT_SOURCE_DIR}/src/hblas/fortran/hblas/hblas1/hdotc.f
    ${PROJECT_SOURCE_DIR}/src/hblas/fortran/hblas/hblas1/hdotu.f
  )
endif()

add_executable(dot ${DOT_SRC})


# Reference BLAS / LAPACK
include(ExternalProject)
find_program(SED_EXEC sed)
ExternalProject_Add(lapack
  PREFIX ${PROJECT_BINARY_DIR}/deps/lapack
  URL "http://www.netlib.org/lapack/lapack-3.5.0.tgz"
  PATCH_COMMAND ${SED_EXEC} -i -e "s/-fltconsistency -fp_port//g" ${PROJECT_BINARY_DIR}/deps/lapack/src/lapack/CMakeLists.txt 
  CMAKE_ARGS -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} 
             -DCMAKE_Fortran_FLAGS='${CMAKE_Fortran_FLAGS}'
             -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/deps
             -DBUILD_TESTING=OFF
)

ExternalProject_Get_Property(lapack install_dir)
add_library(blas STATIC IMPORTED)
set_property(TARGET blas PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/deps/lib/libblas.a)

add_dependencies(blas lapack)
target_link_libraries(gemm ${Boost_LIBRARIES} hblas blas)
target_link_libraries(dot ${Boost_LIBRARIES} hblas)

#find_package(BLAS REQUIRED)
#target_link_libraries(gemm ${Boost_LIBRARIES} haxx_fortran ${BLAS_LIBRARIES})

