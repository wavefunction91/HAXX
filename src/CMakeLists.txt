#
#  This file is a part of HAXX
#
#  Copyright (c) 2017 David Williams-Young
#  All rights reserved.
#
#  See LICENSE.txt 
#

add_subdirectory(hblas)

if( HAXX_ENABLE_BENCHMARK )
  add_subdirectory(benchmark)
endif()
