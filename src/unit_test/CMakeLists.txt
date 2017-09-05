#
#  This file is a part of HAXX
#
#  Copyright (c) 2017 David Williams-Young
#  All rights reserved.
#
#  See LICENSE.txt 
#

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(Boost_USE_STATIC_LIBS ON)
endif()

message( STATUS "HAXX WILL BUILD TEST EXECUTABLES" )
find_package (Boost COMPONENTS unit_test_framework REQUIRED)

add_executable(haxx_unary_op haxx_unary_op.cxx)
add_executable(haxx_binary_op haxx_binary_op.cxx)
add_executable(haxx_algebra haxx_algebra.cxx)
add_executable(haxx_hblas1 haxx_hblas1.cxx)
add_executable(haxx_hblas2 haxx_hblas2.cxx)
add_executable(haxx_hblas3 haxx_hblas3.cxx)
add_executable(haxx_util haxx_util.cxx)

target_link_libraries(haxx_unary_op ${Boost_LIBRARIES} haxx_fortran)
target_link_libraries(haxx_binary_op ${Boost_LIBRARIES} haxx_fortran)
target_link_libraries(haxx_algebra ${Boost_LIBRARIES} haxx_fortran)
target_link_libraries(haxx_hblas1 ${Boost_LIBRARIES} haxx_fortran)
target_link_libraries(haxx_hblas2 ${Boost_LIBRARIES} haxx_fortran)
target_link_libraries(haxx_hblas3 ${Boost_LIBRARIES} haxx_fortran)
target_link_libraries(haxx_util ${Boost_LIBRARIES} haxx_fortran)


add_test(HAXX_UNARY_OP haxx_unary_op --report_level=detailed )  
add_test(HAXX_BINARY_OP haxx_binary_op --report_level=detailed )  
add_test(HAXX_ALGEBRA haxx_algebra --report_level=detailed )  
add_test(HAXX_HBLAS1 haxx_hblas1 --report_level=detailed )  
add_test(HAXX_HBLAS2 haxx_hblas2 --report_level=detailed )  
add_test(HAXX_HBLAS3 haxx_hblas3 --report_level=detailed )  
add_test(HAXX_UTIL haxx_util --report_level=detailed )  
