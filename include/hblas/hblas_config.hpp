/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_CONFIG_HPP
#define __INCLUDED_HBLAS_CONFIG_HPP

#ifndef HAXX_INT
  #define HAXX_INT int32_t
#endif

#define RANK2_INDX(i,j,N) ( (i) + (j)*(N) )

// Hardcoded unrolling parameters
#define HAXX_SWAP_UNROLL 3
#define HAXX_SCAL_UNROLL 5
#define HAXX_COPY_UNROLL 7
#define HAXX_AXPY_UNROLL 4

#include <cassert>

#endif
