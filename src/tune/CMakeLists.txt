#
#  This file is a part of HAXX
#
#  Copyright (c) 2017 David Williams-Young
#  All rights reserved.
#
#  See LICENSE.txt 
#
#set(Boost_USE_STATIC_LIBS ON)



add_executable(gemm_tune gemm.cxx)
target_link_libraries(gemm_tune ${Boost_LIBRARIES} hblas)

