/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#include "haxx.hpp"
#include "hblas/hblas1_def.hpp"

#ifdef __AVX__
#  include <immintrin.h>
#endif

namespace HAXX {

// Optimized Dot Product
  
template <>
quaternion<double> HBLAS_DOTU(const HAXX_INT N, quaternion<double> * const X, 
  const HAXX_INT INCX, quaternion<double> * const Y, const HAXX_INT INCY) {

  quaternion<double> res(0.);

  quaternion<double> * locX = X, *locY = Y;

#ifdef __AVX__

  // Load quaternions
  __m256d x1;
  __m256d x2;
  __m256d x3;
  __m256d x4;

  __m256d y1;
  __m256d y2;
  __m256d y3;
  __m256d y4;
  

  // Scratch space
  __m256d t1, t2, t3, t4;

  // Result buffer
  __m256d r1 = _mm256_setzero_pd();
  __m256d r2 = _mm256_setzero_pd();
  __m256d r3 = _mm256_setzero_pd();
  __m256d r4 = _mm256_setzero_pd();
  
#endif

  const HAXX_INT BatchSize = 4;
  const HAXX_INT NBatch    = N / BatchSize;
  const HAXX_INT nLeft = N % BatchSize;

  HAXX_INT i;


  const bool XIsAligned = ((unsigned long)X & 31) == 0;
  const bool YIsAligned = ((unsigned long)Y & 31) == 0;

  const bool isAligned = XIsAligned and YIsAligned;


  for( i = 0; i < NBatch; i++ ) {
  
#ifdef __AVX__
            
    if( isAligned ) {

      x1 = _mm256_load_pd( reinterpret_cast<double*>(locX) );
      x2 = _mm256_load_pd( reinterpret_cast<double*>(locX + INCX  ) );
      x3 = _mm256_load_pd( reinterpret_cast<double*>(locX + 2*INCX) );
      x4 = _mm256_load_pd( reinterpret_cast<double*>(locX + 3*INCX) );

      y1 = _mm256_load_pd( reinterpret_cast<double*>(locY) );
      y2 = _mm256_load_pd( reinterpret_cast<double*>(locY + INCY) );
      y3 = _mm256_load_pd( reinterpret_cast<double*>(locY + 2*INCY) );
      y4 = _mm256_load_pd( reinterpret_cast<double*>(locY + 3*INCY) );

    } else {

      x1 = _mm256_loadu_pd( reinterpret_cast<double*>(locX) );
      x2 = _mm256_loadu_pd( reinterpret_cast<double*>(locX + INCX  ) );
      x3 = _mm256_loadu_pd( reinterpret_cast<double*>(locX + 2*INCX) );
      x4 = _mm256_loadu_pd( reinterpret_cast<double*>(locX + 3*INCX) );

      y1 = _mm256_loadu_pd( reinterpret_cast<double*>(locY) );
      y2 = _mm256_loadu_pd( reinterpret_cast<double*>(locY + INCY) );
      y3 = _mm256_loadu_pd( reinterpret_cast<double*>(locY + 2*INCY) );
      y4 = _mm256_loadu_pd( reinterpret_cast<double*>(locY + 3*INCY) );

    }


    // Transpose X batch
    t1 = _mm256_shuffle_pd(x1, x2, 0x0);
    t3 = _mm256_shuffle_pd(x1, x2, 0xF);
    t2 = _mm256_shuffle_pd(x3, x4, 0x0);
    t4 = _mm256_shuffle_pd(x3, x4, 0xF);

    x1 = _mm256_permute2f128_pd(t1, t2, 0x20);
    x2 = _mm256_permute2f128_pd(t3, t4, 0x20);
    x3 = _mm256_permute2f128_pd(t1, t2, 0x31);
    x4 = _mm256_permute2f128_pd(t3, t4, 0x31);

    // Transpose Y batch
    t1 = _mm256_shuffle_pd(y1, y2, 0x0);
    t3 = _mm256_shuffle_pd(y1, y2, 0xF);
    t2 = _mm256_shuffle_pd(y3, y4, 0x0);
    t4 = _mm256_shuffle_pd(y3, y4, 0xF);

    y1 = _mm256_permute2f128_pd(t1, t2, 0x20);
    y2 = _mm256_permute2f128_pd(t3, t4, 0x20);
    y3 = _mm256_permute2f128_pd(t1, t2, 0x31);
    y4 = _mm256_permute2f128_pd(t3, t4, 0x31);


    // SCALAR parts of product
    // R(S) = X(S)*Y(S) - X(I)*Y(I) - X(J)*Y(J) - X(K)*Y(K)

    // t1 = SCALAR * SCALAR
    t1 = _mm256_mul_pd(x1,y1);


    // t1 = t1 - I * I
    t2 = _mm256_mul_pd(x2,y2);
    t1 = _mm256_sub_pd(t1,t2);

    
    // t1 = t1 - J * J
    t2 = _mm256_mul_pd(x3,y3);
    t1 = _mm256_sub_pd(t1,t2);

    
    // r1 += t1 - K * K
    t2 = _mm256_mul_pd(x4,y4);
    r1 = _mm256_add_pd(r1,t1);
    r1 = _mm256_sub_pd(r1,t2);




    // I parts of product

    // t1 = SCALAR * I
    t1 = _mm256_mul_pd(x1,y2);
   
    // t1 = t1 + I * SCALAR      
    t2 = _mm256_mul_pd(x2,y1);
    t1 = _mm256_add_pd(t1,t2);

    // t1 = t1 + J * K
    t2 = _mm256_mul_pd(x3,y4);
    t1 = _mm256_add_pd(t1,t2);

    // r2 = t1 - K * J
    t2 = _mm256_mul_pd(x4,y3);
    r2 = _mm256_add_pd(r2,t1);
    r2 = _mm256_sub_pd(r2,t2);


    // J parts of product

    // t1 = SCALAR * J
    t1 = _mm256_mul_pd(x1,y3);
   
    // t1 = t1 - I * K      
    t2 = _mm256_mul_pd(x2,y4);
    t1 = _mm256_sub_pd(t1,t2);

    // t1 = t1 + J * SCALAR
    t2 = _mm256_mul_pd(x3,y1);
    t1 = _mm256_add_pd(t1,t2);

    // t1 = t1 + K * I
    t2 = _mm256_mul_pd(x4,y2);
    r3 = _mm256_add_pd(r3,t1);
    r3 = _mm256_add_pd(r3,t2);


    // K parts of product

    // t1 = SCALAR * K
    t1 = _mm256_mul_pd(x1,y4);
   
    // t1 = t1 + I * J      
    t2 = _mm256_mul_pd(x2,y3);
    t1 = _mm256_add_pd(t1,t2);

    // t1 = t1 - J * I
    t2 = _mm256_mul_pd(x3,y2);
    t1 = _mm256_sub_pd(t1,t2);

    // t1 = t1 + K * SCALAR
    t2 = _mm256_mul_pd(x4,y1);
    r4 = _mm256_add_pd(r4,t1);
    r4 = _mm256_add_pd(r4,t2);


#endif


    locX += BatchSize * INCX;
    locY += BatchSize * INCY;
  
  }


  double * resPtr = (double*)(&res);

#ifdef __AVX__

  // Transpose R batch
  t1 = _mm256_shuffle_pd(r1, r2, 0x0);
  t3 = _mm256_shuffle_pd(r1, r2, 0xF);
  t2 = _mm256_shuffle_pd(r3, r4, 0x0);
  t4 = _mm256_shuffle_pd(r3, r4, 0xF);

  r1 = _mm256_permute2f128_pd(t1, t2, 0x20);
  r2 = _mm256_permute2f128_pd(t3, t4, 0x20);
  r3 = _mm256_permute2f128_pd(t1, t2, 0x31);
  r4 = _mm256_permute2f128_pd(t3, t4, 0x31);

  r1 = _mm256_add_pd(r1,r2);
  r1 = _mm256_add_pd(r1,r3);
  r1 = _mm256_add_pd(r1,r4);


  _mm256_storeu_pd(resPtr,r1);   

#endif

  for(i = 0; i < nLeft; i++)
    res = res + locX[i*INCX] * locY[i*INCY];
    
  return res;

}; // HBLAS_DOTU

}; // namespace HAXX
