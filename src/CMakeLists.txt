#
#  This file is a part of HAXX
#
#  Copyright (c) 2017 David Williams-Young
#  All rights reserved.
#
#  See LICENSE.txt 
#
#add_subdirectory(scalar_ops)

add_subdirectory(fortran_bindings)

if( HAXX_ENABLE_TESTING )
  add_subdirectory(unit_test)
endif()

if( HAXX_ENABLE_BENCHMARK )
  add_subdirectory(benchmark)
endif()
