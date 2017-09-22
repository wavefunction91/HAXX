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

#include <iostream>

#include <util/simd.hpp>


// Determine type of scaling parameter ALPHA
#if ALPHAF == DOUBLE
  #define _ALPHAF double
#elif ALPHAF == DCOMPLEX
  #define _ALPHAF std::complex<double>
#elif ALPHAF == DQUATERNION
  #define _ALPHAF quaternion<double>
#else
  #error AXPY/SCAL Only Supports 64-bit floats
#endif


namespace HAXX {

// Optimized AXPY and SCAL operations

template <>
#ifdef _AXPY

void FNAME(const char SIDE, const HAXX_INT N, const _ALPHAF ALPHA, 
  quaternion<double> * const X, const HAXX_INT INCX, 
  quaternion<double> * const Y, const HAXX_INT INCY)

#elif defined(_SCAL)

void FNAME(const char SIDE, const HAXX_INT N, const _ALPHAF ALPHA, 
  quaternion<double> * const X, const HAXX_INT INCX)

#endif
{

  // Local vector pointers for increment
  quaternion<double> * locX = X;

#ifdef _AXPY
  quaternion<double> * locY = Y;
#endif

#if defined(__AVX__) || defined(__AVX2__)

  // Load quaternions
  __m256d x1; // Load space for an element of X

  #ifdef _AXPY
  __m256d y1; // Load space for an element of Y
  #endif
  
#endif

  const HAXX_INT BatchSize = 1;
  const HAXX_INT NBatch    = N / BatchSize;
  const HAXX_INT nLeft = N % BatchSize;

  HAXX_INT i;

  // Determine side
  const bool MulLeft = (SIDE == 'L');

  // Determine alignment
  const bool XIsAligned = IS_ALIGNED(X,REQ_ALIGN);

#ifdef _AXPY
  const bool YIsAligned = IS_ALIGNED(Y,REQ_ALIGN);
  const bool isAligned  = XIsAligned and YIsAligned;
#else
  const bool isAligned  = XIsAligned;
#endif

  // Load scaling factor
#if defined(__AVX__) || defined(__AVX2__)


  #if ALPHAF == DOUBLE

  // alpha = (ALPHA, ALPHA, ALPHA, ALPHA)
  const __m256d alpha = _mm256_broadcast_sd(&ALPHA);

  #elif ALPHAF == DCOMPLEX


  // Quaternion product with a complex number may be thought of as
  // two complex products via tha Cayley-Dickson construction of the
  // quaternions
  //
  // A * (B + Cj) = AB + ACj
  // (A + Bj) * C = AC + B*CONJ(C) j


  auto ALPHA_C = std::conj(ALPHA);

  double *ALPHA_PTR = const_cast<double*>(
                        reinterpret_cast<const double *>(&ALPHA));

  __m128d alphaC      = _mm_loadu_pd(ALPHA_PTR);


  // Precompute the second vector for complex multiplication
  __m128d alphaC_conj = 
    MulLeft ? _mm_loadu_pd(ALPHA_PTR) : 
              _mm_loadu_pd(reinterpret_cast<double*>(&ALPHA_C));



  // SIDE == 'L' -> 
  //   alpha      = (ALPHA_R,  ALPHA_I, ALPHA_R,  ALPHA_I)
  //   alpha_conj = (ALPHA_I, -ALPHA_R, ALPHA_I, -ALPHA_R)
  // SIDE == 'R' -> 
  //   alpha      = (ALPHA_R,  ALPHA_I,  ALPHA_R, -ALPHA_I)
  //   alpha_conj = (ALPHA_I, -ALPHA_R, -ALPHA_I, -ALPHA_R)
  const __m256d alpha = D256_FROM_D128(_mm_castpd_ps(alphaC),
                                    _mm_castpd_ps(alphaC_conj));

  // (0 -1 0 -1) mask
  const __m256i maskConj = _mm256_set_epi64x(
                             0x8000000000000000, 0,
                             0x8000000000000000, 0 );

  const __m256d alpha_conj = _mm256_xor_pd(
                            _mm256_permute_pd(alpha, 0x5),
                            _mm256_castsi256_pd(maskConj)
                          );

  #elif ALPHAF == DQUATERNION

  // Load ALPHA as is
  const __m256d alpha = LOADD_UNALIGNED_AS(double,&ALPHA);

  #endif

#endif

  for( i = 0; i < NBatch; i++ ) {

#if defined(__AVX__) || defined (__AVX2__)

    // Load X and possibly Y
    if( isAligned ) {

      x1 = LOADD_ALIGNED_AS(double,locX);
    #ifdef _AXPY
      y1 = LOADD_ALIGNED_AS(double,locY);
    #endif

    } else {

      x1 = LOADD_UNALIGNED_AS(double,locX);
    #ifdef _AXPY
      y1 = LOADD_UNALIGNED_AS(double,locY);
    #endif

    }

#ifdef _SCAL

  // X = ALPHA * X or X * ALPHA
  #if   ALPHAF == DOUBLE
    x1 = MULD(alpha,x1);
  #elif ALPHAF == DCOMPLEX
    __m256d p1 = MULD(x1,alpha);
    __m256d p2 = MULD(x1,alpha_conj);
    x1      = _mm256_hsub_pd(p1,p2);
  #elif ALPHAF == DQUATERNION
    if( MulLeft ) x1 = MULDQ_NN(alpha,x1);
    else          x1 = MULDQ_NN(x1,alpha);
  #endif

  // Store X in place
  if( isAligned ) STORED_ALIGNED_AS(  double,locX,x1);
  else            STORED_UNALIGNED_AS(double,locX,x1);        

#elif defined(_AXPY)

  // Y = Y + ALPHA * X or Y + X * ALPHA
  #if   ALPHAF == DOUBLE
    y1 = FMAD(alpha,x1,y1);
  #elif ALPHAF == DCOMPLEX
    __m256d p1 = MULD(x1,alpha);
    __m256d p2 = MULD(x1,alpha_conj);
    y1      = ADDD(y1,_mm256_hsub_pd(p1,p2));
  #elif ALPHAF == DQUATERNION
    if( MulLeft ) y1 = ADDD(y1,MULDQ_NN(alpha,x1));
    else          y1 = ADDD(y1,MULDQ_NN(x1,alpha));
  #endif

  // Store Y in place
  if( isAligned ) STORED_ALIGNED_AS(  double,locY,y1);
  else            STORED_UNALIGNED_AS(double,locY,y1);        

  #endif

#endif

    // Increment vectors
    locX += BatchSize * INCX;
  #ifdef _AXPY
    locY += BatchSize * INCY;
  #endif

  }

  // FIXME: this only works for a batch size of 1, need a
  // cleanup loop for AVX-512
  //
};

}; // namespace HAXX
