/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_SIMD_MISC_HPP__
#define __INCLUDED_SIMD_MISC_HPP__


// Misc macros


// Alignment checking
#define IS_ALIGNED(X,B) ( ((unsigned long)(X) & (B-1)) == 0 )


// Load operations

// Load aligned with cast
#define LOADD_ALIGNED_AS(T,X) \
  LOADD_ALIGNED(const_cast<T*>(reinterpret_cast<const T*>(X)))

// Load unaligned with cast
#define LOADD_UNALIGNED_AS(T,X) \
  LOADD_UNALIGNED(const_cast<T*>(reinterpret_cast<const T*>(X)))

// Store aligned with cast
#define STORED_ALIGNED_AS(T,X,V) \
  STORED_ALIGNED(reinterpret_cast<T*>(X),V)

// Load unaligned with cast
#define STORED_UNALIGNED_AS(T,X,V) \
  STORED_UNALIGNED(reinterpret_cast<T*>(X),V)




#ifdef __AVX__

  // Transpose 4x4 registers (with scratch space)
  #define _MM_TRANSPOSE_4x4_PD(w,x,y,z,t1,t2,t3,t4) \
    t1 = _mm256_shuffle_pd(w, x, 0x0);\
    t3 = _mm256_shuffle_pd(w, x, 0xF);\
    t2 = _mm256_shuffle_pd(y, z, 0x0);\
    t4 = _mm256_shuffle_pd(y, z, 0xF);\
    \
    w = _mm256_permute2f128_pd(t1, t2, 0x20);\
    x = _mm256_permute2f128_pd(t3, t4, 0x20);\
    y = _mm256_permute2f128_pd(t1, t2, 0x31);\
    z = _mm256_permute2f128_pd(t3, t4, 0x31);
    
  
  #endif

#endif
