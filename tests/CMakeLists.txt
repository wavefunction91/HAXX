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

# HAXX UTs
add_executable(haxx_ut ut.cxx haxx_algebra.cxx haxx_binary_op.cxx 
  haxx_unary_op.cxx)
target_compile_definitions(haxx_ut PUBLIC BOOST_TEST_MODULE=HAXX_FUNCTIONALITY)
target_link_libraries(haxx_ut ${Boost_LIBRARIES})

add_test(HAXX_ALGEBRA   haxx_ut --report_level=detailed --run_test=HAXX_ALGEBRA  )
add_test(HAXX_BINARY_OP haxx_ut --report_level=detailed --run_test=HAXX_BINARY_OP)
add_test(HAXX_UNARY_OP  haxx_ut --report_level=detailed --run_test=HAXX_UNARY_OP )

add_executable(hblas_ut ut.cxx hblas1.cxx hblas2.cxx hblas3.cxx hblas_util.cxx)
target_compile_definitions(hblas_ut PUBLIC BOOST_TEST_MODULE=HBLAS_FUNCTIONALITY)
target_link_libraries(hblas_ut ${Boost_LIBRARIES} hblas)

add_test(HBLAS_HBLAS1 hblas_ut --report_level=detailed --run_test=HBLAS1)
add_test(HBLAS_HBLAS2 hblas_ut --report_level=detailed --run_test=HBLAS2)
add_test(HBLAS_HBLAS3 hblas_ut --report_level=detailed --run_test=HBLAS3)
add_test(HBLAS_UTIL   hblas_ut --report_level=detailed --run_test=UTIL  )
