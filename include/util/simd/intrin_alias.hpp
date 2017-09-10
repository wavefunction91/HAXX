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
#ifdef __AVX__

  // 64-bit float SIMD vector
  #define VECD __m256d


  // 64-bit float Zero vector
  #define ZEROD _mm256_setzero_pd()



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

#endif




// FMA / FMS operations (D = A*B +- C)

// AVX has no FMA, macro emulates with 2 instructions
#ifdef __AVX__

#define FMAD(A,B,C) ADDD(MULD(A,B),C)
#define FMSD(A,B,C) SUBD(MULD(A,B),C)

#endif




// Macros for FSM (D = A +- B*C)
// XXX: This is not a single FMA type instruction

#ifdef __AVX__
#  define FSMD(A,B,C) SUBD(A,MULD(B,C))
#endif




#endif
