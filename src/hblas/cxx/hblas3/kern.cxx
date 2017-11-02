/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#include "haxx.hpp"
#include "hblas/hblas1.hpp"
#include "hblas/hblas3.hpp"

#include "util/simd.hpp"
#include "util/macro.hpp"

#include "hblas/config/types.hpp"
#include "hblas/config/hblas3/gemm.hpp"


namespace HAXX {

template <typename T, typename U, typename V>
 void Kern(HAXX_INT M, HAXX_INT N, HAXX_INT K, T* __restrict__ A, U* __restrict__ B, V* __restrict__ C, 
  HAXX_INT LDC);

template<>
 void Kern(HAXX_INT M, HAXX_INT N, HAXX_INT K, 
  quaternion<double>* __restrict__ A, quaternion<double>* __restrict__ B, 
  quaternion<double>* __restrict__ C, HAXX_INT LDC) {

  __m256d t1,t2,t3,t4;

  const bool M2 = (M == 2);
  const bool N2 = (N == 2);

  // Load C
  __m256d c00 = LOAD_256D_UNALIGNED_AS(double,C      );

  __m256d c10 = M2 ? LOAD_256D_UNALIGNED_AS(double,C+1    ) :
                           _mm256_setzero_pd();

  __m256d c01 = N2 ? LOAD_256D_UNALIGNED_AS(double,C+LDC  ) :
                           _mm256_setzero_pd();

  __m256d c11 = ( M2 and N2 ) ?  LOAD_256D_UNALIGNED_AS(double,C+LDC+1) :
                                 _mm256_setzero_pd();

  _MM_TRANSPOSE_4x4_PD(c00,c01,c10,c11,t1,t2,t3,t4);


  quaternion<double> *locB = B, *locA = A;

  HAXX_INT k = K;

  if( k > 0 )
  do {

    // Load A 
    __m256d a00 = LOAD_256D_ALIGNED_AS(double,locA);
    __m256d a10 = LOAD_256D_ALIGNED_AS(double,locA+1);
    locA += 2;
    

    // Load B

#ifndef _FACTOR_TRANSPOSE_INTO_B_PACK
    __m256d b00 = LOAD_256D_ALIGNED_AS(double,locB);
    __m256d b10 = LOAD_256D_ALIGNED_AS(double,locB+1);
#else
    double *BasDouble = reinterpret_cast<double*>(locB);
    __m128d b00lo = LOAD_128D_ALIGNED(BasDouble);
    __m128d b00hi = LOAD_128D_ALIGNED(BasDouble+2);
    __m128d b10lo = LOAD_128D_ALIGNED(BasDouble+4);
    __m128d b10hi = LOAD_128D_ALIGNED(BasDouble+6);
#endif
    locB += 2;


#ifdef _FACTOR_TRANSPOSE_INTO_A_PACK

    __m256d a_IIII = _mm256_permute_pd(a00,0xF);
    a00            = _mm256_permute_pd(a00,0x0); // SSSS
  
    __m256d a_KKKK = _mm256_permute_pd(a10,0xF);
    a10            = _mm256_permute_pd(a10,0x0); // SSSS

    __m256d &a00c = a_IIII;
    __m256d &a10c = a_KKKK;

#else

    __m256d a00c = a00;
    __m256d a10c = a10;
    _MM_TRANSPOSE_4x4_PD(a00,a00c,a10,a10c,t1,t2,t3,t4);

#endif

#ifdef _FACTOR_TRANSPOSE_INTO_B_PACK

    __m256d bSSSS = SET_256D_FROM_128D(b00lo,b00lo);
    __m256d bIIII = SET_256D_FROM_128D(b10lo,b10lo);
    __m256d bJJJJ = SET_256D_FROM_128D(b00hi,b00hi);
    __m256d bKKKK = SET_256D_FROM_128D(b10hi,b10hi);

    __m256d &b00  = bSSSS;
    __m256d &b10  = bIIII;
    __m256d &b00c = bJJJJ;
    __m256d &b10c = bKKKK;

#else

    __m256d b00c = b00;
    __m256d b10c = b10;
    _MM_TRANSPOSE_4x4_PD(b00,b10,b00c,b10c,t1,t2,t3,t4);

#endif

    INC_MULD4Q_NN(a00,a00c,a10,a10c,b00,b10,b00c,b10c,c00,c01,c10,c11);
  

    k--;
  } while( k > 0 );

  _MM_TRANSPOSE_4x4_PD(c00,c01,c10,c11,t1,t2,t3,t4);

  STORE_256D_UNALIGNED_AS(double,C      ,c00);
  if( M2 ) STORE_256D_UNALIGNED_AS(double,C+1    ,c10);
  if( N2 ) STORE_256D_UNALIGNED_AS(double,C+LDC  ,c01);
  if( M2 and N2 ) STORE_256D_UNALIGNED_AS(double,C+LDC+1,c11);
  

}
};
