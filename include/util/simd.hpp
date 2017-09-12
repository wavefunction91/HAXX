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
  
#include <immintrin.h>

// Compile time type constants
#define DOUBLE      0
#define DCOMPLEX    1
#define DQUATERNION 2

// Required boundary alignment for aligned data
#if defined(__AVX__) || defined(__AVX2__)
  #define REQ_ALIGN 32
#endif



#include "simd/intrin_alias.hpp"
#include "simd/misc.hpp"
#include "simd/qop.hpp"


#endif
