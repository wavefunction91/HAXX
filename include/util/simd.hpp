/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_UTIL_SIMD_HPP__
#define __INCLUDED_UTIL_SIMD_HPP__

// SIMD intrinsics header
  
#ifdef __AVX__
  #include <immintrin.h>
#endif

// Compile time type constants
#define DOUBLE      0
#define DCOMPLEX    1
#define DQUATERNION 2

// Required boundary alignment for aligned data
#ifdef __AVX__
  #define REQ_ALIGN 32
#endif



#include "simd/intrin_alias.hpp"
#include "simd/misc.hpp"
#include "simd/qop.hpp"


#endif
