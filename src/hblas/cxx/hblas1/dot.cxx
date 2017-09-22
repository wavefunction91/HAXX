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

#include "util/simd.hpp"

namespace HAXX {

// Optimized Dot Product
  
template <>
quaternion<double> FNAME(const HAXX_INT N, quaternion<double> * const X, 
  const HAXX_INT INCX, quaternion<double> * const Y, const HAXX_INT INCY) {


  // Result quaternion
  quaternion<double> res(0.);

  // Local pointers for increment
  quaternion<double> * locX = X, *locY = Y;

#if defined(__AVX__) || defined(__AVX2__)

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

  // Result buffer (zeroed out)
  __m256d r1 = ZEROD;
  __m256d r2 = ZEROD;
  __m256d r3 = ZEROD;
  __m256d r4 = ZEROD;
  
#endif

  // Batch size information
  const HAXX_INT BatchSize = 4;
  const HAXX_INT NBatch    = N / BatchSize;
  const HAXX_INT nLeft = N % BatchSize;

  HAXX_INT i;


  // Determine alignment
  const bool XIsAligned = IS_ALIGNED(X,REQ_ALIGN);
  const bool YIsAligned = IS_ALIGNED(Y,REQ_ALIGN);

  const bool isAligned = XIsAligned and YIsAligned;

  for( i = 0; i < NBatch; i++ ) {
  
#if defined(__AVX__) || defined(__AVX2__)
            
    // Load 4 X and Y elements
    if( isAligned ) {

      x1 = LOADD_ALIGNED_AS(double,locX         );
      x2 = LOADD_ALIGNED_AS(double,locX +   INCX);
      x3 = LOADD_ALIGNED_AS(double,locX + 2*INCX);
      x4 = LOADD_ALIGNED_AS(double,locX + 3*INCX);

      y1 = LOADD_ALIGNED_AS(double,locY         );
      y2 = LOADD_ALIGNED_AS(double,locY +   INCY);
      y3 = LOADD_ALIGNED_AS(double,locY + 2*INCY);
      y4 = LOADD_ALIGNED_AS(double,locY + 3*INCY);

    } else {

      x1 = LOADD_UNALIGNED_AS(double,locX         );
      x2 = LOADD_UNALIGNED_AS(double,locX +   INCX);
      x3 = LOADD_UNALIGNED_AS(double,locX + 2*INCX);
      x4 = LOADD_UNALIGNED_AS(double,locX + 3*INCX);

      y1 = LOADD_UNALIGNED_AS(double,locY         );
      y2 = LOADD_UNALIGNED_AS(double,locY +   INCY);
      y3 = LOADD_UNALIGNED_AS(double,locY + 2*INCY);
      y4 = LOADD_UNALIGNED_AS(double,locY + 3*INCY);

    }


    // Transpose X batch
    _MM_TRANSPOSE_4x4_PD(x1,x2,x3,x4,t1,t2,t3,t4);

    // Transpose Y batch
    _MM_TRANSPOSE_4x4_PD(y1,y2,y3,y4,t1,t2,t3,t4);

    // R += CONJ(X) * Y
    #ifdef _CONJ
      INC_MULD4Q_CN(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);

    // R += X * Y
    #else
      INC_MULD4Q_NN(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);

    #endif



#endif


    // Increment pointers
    locX += BatchSize * INCX;
    locY += BatchSize * INCY;
  
  }


#if defined(__AVX__) || defined(__AVX2__)

  // Transpose R batch to proper quaternion storage
  _MM_TRANSPOSE_4x4_PD(r1,r2,r3,r4,t1,t2,t3,t4);

  // Add up the four quaternions
  // r1 = r1 + r2 + r3 + r4
  r1 = ADDD(r1,r2);
  r1 = ADDD(r1,r3);
  r1 = ADDD(r1,r4);

#endif

  // Cleanup loop, multiply one quaternion at a time
  for(i = 0; i < nLeft; i++) {

    // Load a single element of X and Y
    if( isAligned ) {
      x1 = LOADD_ALIGNED_AS(double,locX + i*INCX);
      y1 = LOADD_ALIGNED_AS(double,locY + i*INCY);
    } else {
      x1 = LOADD_UNALIGNED_AS(double,locX + i*INCX);
      y1 = LOADD_UNALIGNED_AS(double,locY + i*INCY);
    }

    // r1 += CONJ(X) * Y
  #ifdef _CONJ
    r1 = ADDD(r1,MULDQ_CN(x1,y1));

    // r1 += X * Y
  #else
    r1 = ADDD(r1,MULDQ_NN(x1,y1));

  #endif

  } 

  // Store the result in a possibly unaligned
  // storage
  STORED_UNALIGNED_AS(double,(&res),r1);   

    
  // return
  return res;

}; // HBLAS_DOT?

}; // namespace HAXX
