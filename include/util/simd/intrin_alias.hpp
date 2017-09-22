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


// 256-bit vectors

// Load operations
#define LOAD_256D_ALIGNED(X)    _mm256_load_pd(X)
#define LOAD_256D_UNALIGNED(X)  _mm256_loadu_pd(X)

// Load operations with proper cast
#define LOAD_256D_ALIGNED_AS(T,X)\
  LOAD_256D_ALIGNED(const_cast<T*>(reinterpret_cast<const T*>(X)))
#define LOAD_256D_UNALIGNED_AS(T,X)\
  LOAD_256D_UNALIGNED(const_cast<T*>(reinterpret_cast<const T*>(X)))


// Store operations
#define STORE_256D_ALIGNED(X,V)    _mm256_store_pd(X,V)
#define STORE_256D_UNALIGNED(X,V)  _mm256_storeu_pd(X,V)

// Store operations with proper cast
#define STORE_256D_ALIGNED_AS(T,X,V)\
  STORE_256D_ALIGNED(reinterpret_cast<T*>(X),V)
#define STORE_256D_UNALIGNED_AS(T,X,V)\
  STORE_256D_UNALIGNED(reinterpret_cast<T*>(X),V)


// Getting / Assembling operations
#define GET_LO_128D_256D(X) _mm256_castpd256_pd128(X)
#define GET_HI_128D_256D(X) _mm256_extractf128_pd((X),1)

#define SET_256D_FROM_128D(X,Y) \
  _mm256_castps_pd(\
    _mm256_insertf128_ps(\
      _mm256_castps128_ps256(_mm_castpd_ps(X)),\
      _mm_castpd_ps(Y),1)\
  );













// FMA / FMS operations (D = +-(A*B) +- C)

// Emulate FMA with two instructions
#ifndef __FMA__

  // 256-bit vectors
  #define FMA_256D(A,B,C)  _mm256_add_pd(_mm256_mul_pd(A,B),C)
  #define FMS_256D(A,B,C)  _mm256_sub_pd(_mm256_mul_pd(A,B),C)
  #define FNMA_256D(A,B,C) _mm256_sub_pd(C,_mm256_mul_pd(A,B))

#else

  // 256-bit vectors
  #define FMA_256D(A,B,C)  _mm256_fmadd_pd(A,B,C)
  #define FMS_256D(A,B,C)  _mm256_fmsub_pd(A,B,C)
  #define FNMA_256D(A,B,C) _mm256_fnmadd_pd(A,B,C)

#endif












#endif
