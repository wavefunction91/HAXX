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

// Optimized Copy

template <>
void FNAME(const HAXX_INT N, quaternion<double> * const X, 
  const HAXX_INT INCX, quaternion<double> * const Y, const HAXX_INT INCY) {

  // Local pointers for increment
  quaternion<double> * locX = X, *locY = Y;

#if defined(__AVX__) || defined(__AVX2__)

  __m256d x1;

#ifdef _SWAP

  __m256d y1;

#endif

#endif

  // Batch size information
  const HAXX_INT BatchSize = 1;
  const HAXX_INT NBatch    = N / BatchSize;
  const HAXX_INT nLeft = N % BatchSize;

  HAXX_INT i;

  assert(BatchSize == 1); // Need to address increment for AVX512

  // Determine alignment
  const bool XIsAligned = IS_ALIGNED(X,REQ_ALIGN);
  const bool YIsAligned = IS_ALIGNED(Y,REQ_ALIGN);

  const bool isAligned = XIsAligned and YIsAligned;


  if( isAligned )
    for( i = 0; i < NBatch; i++ ) {

  #if defined(__AVX__) || defined(__AVX2__)

      x1 = LOADD_ALIGNED_AS(double,locX);
    #ifdef _SWAP
      y1 = LOADD_ALIGNED_AS(double,locY);
    #endif

      STORED_ALIGNED_AS(double,locY,x1);
    #ifdef _SWAP
      STORED_ALIGNED_AS(double,locX,y1);
    #endif

  #endif

      locX += INCX;
      locY += INCY;

    }
  else
    for( i = 0; i < NBatch; i++ ) {
  
  #if defined(__AVX__) || defined(__AVX2__)
  
      x1 = LOADD_UNALIGNED_AS(double,locX);
    #ifdef _SWAP
      y1 = LOADD_UNALIGNED_AS(double,locY);
    #endif


      STORED_UNALIGNED_AS(double,locY,x1);
    #ifdef _SWAP
      STORED_UNALIGNED_AS(double,locX,y1);
    #endif
  
  #endif
  
      locX += INCX;
      locY += INCY;
  
    }

}; // HBLAS_COPY / HBLAS_SWAP

}; // namspace HAXX
