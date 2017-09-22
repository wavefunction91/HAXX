/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_SIMD_INTRIN_ALIAS_HPP__
#define __INCLUDED_SIMD_INTRIN_ALIAS_HPP__


// Alias SIMD intrinsics

// 256-bit vector length
#if defined(__AVX__) || defined(__AVX2__)

  // Load/Store operations
  
  // Load Aligned
  #define LOADD_ALIGNED(X) _mm256_load_pd(X)
  
  // Load Unaligned
  #define LOADD_UNALIGNED(X) _mm256_loadu_pd(X)
  
  // Load Aligned
  #define STORED_ALIGNED(X,V) _mm256_store_pd(X,V)
  
  // Load Unaligned
  #define STORED_UNALIGNED(X,V) _mm256_storeu_pd(X,V)






  // Arithmetic
  
  // Add 64-bit floats
  #define ADDD(X,Y) _mm256_add_pd(X,Y)
  
  // Subtract 64-bit floats
  #define SUBD(X,Y) _mm256_sub_pd(X,Y)
  
  // Multiply 64-bit floats
  #define MULD(X,Y) _mm256_mul_pd(X,Y)



  // Get the low and high 128-bits of 256-bit register
  #define GET_LOD_256(X) _mm256_castpd256_pd128(X)
  #define GET_HID_256(X) _mm256_extractf128_pd((X),1)


  // Construct a 256-bit pack from 2 128-bit packs
  #define D256_FROM_D128(X,Y) \
    _mm256_castps_pd(_mm256_insertf128_ps(_mm256_castps128_ps256(X),(Y),1));


#endif




// FMA / FMS operations (D = +-(A*B) +- C)

// Emulate FMA with two instructions
#ifndef __FMA__

  #define FMAD(A,B,C)  ADDD(MULD(A,B),C)
  #define FMSD(A,B,C)  SUBD(MULD(A,B),C)
  #define FNMAD(A,B,C) SUBD(C,MULD(A,B))

#else

  #define FMAD(A,B,C)  _mm256_fmadd_pd(A,B,C)
  #define FMSD(A,B,C)  _mm256_fmsub_pd(A,B,C)
  #define FNMAD(A,B,C) _mm256_fnmadd_pd(A,B,C)

#endif




#endif
