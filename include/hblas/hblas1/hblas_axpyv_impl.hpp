/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#pragma once

#include "hblas/hblas1_def.hpp"

namespace HAXX {

/**
 *  Written by DBWY (4/2017)
 *
 *  Based on the BLAS implementaion of DCOPY by Jack Dongarra
 *    http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
 *
 *
 *  \f$ y \in \mathbb{H} \qquad x,\alpha \in \mathbb{R},\mathbb{C},\mathbb{H} \f$
 *
 *  SIDE == 'L'
 *
 *  \f$ y_i = \alpha x_i  + y_i\f$
 *
 *  SIDE == 'R'
 *
 *  \f$ y_i = x_i \alpha + y_i \f$
 */
template <typename _F, typename _XF, typename _AlphaF>
void HBLAS_AXPYV(const char SIDE, const HAXX_INT N, const _AlphaF ALPHA, 
  _XF * const X, const HAXX_INT INCX, quaternion<_F> * const Y, 
  const HAXX_INT INCY) {

  
  if( N <= 0 ) return;
  if( ALPHA == _AlphaF(0.) ) return;

  // FIXME: See further comments on negative stride
  assert(INCX > 0);
  assert(INCY > 0);

  HAXX_INT i;

  // FIXME: Should write a specialization for real ALPHA where side
  //   doesnt matter
  bool isR = SIDE == 'R';
  bool isL = not isR;

  if( INCX == 1 and INCY == 1 ) {
    HAXX_INT m = N % HAXX_AXPY_UNROLL;

    if( m != 0) {
      if( isL ) for( i = 0; i < m; ++i ) Y[i] += ALPHA * X[i];
      else      for( i = 0; i < m; ++i ) Y[i] += X[i] * ALPHA;

      // XXX: DAXPY has this outside of the if-check? Unline COPY and SCAL
      if( N < HAXX_AXPY_UNROLL ) return;
    }

    // FIXME: This assumes HAXX_AXPY_UNROLL = 4
    if( isL )
      for( i = m; i < N; i += HAXX_AXPY_UNROLL ) {
        Y[i]   += ALPHA * X[i];
        Y[i+1] += ALPHA * X[i+1];
        Y[i+2] += ALPHA * X[i+2];
        Y[i+3] += ALPHA * X[i+3];
      }
    else
      for( i = m; i < N; i += HAXX_AXPY_UNROLL ) {
        Y[i]   += X[i]   * ALPHA; 
        Y[i+1] += X[i+1] * ALPHA; 
        Y[i+2] += X[i+2] * ALPHA; 
        Y[i+3] += X[i+3] * ALPHA; 
      }
  } else {

    HAXX_INT ix(0), iy(0);
    // FIXME: the original _AXPY function has code here to handle
    //   negative increments. Unsure on what that accomplishes

    if( isL ) 
      for( i = 0; i < N; ++i, ix += INCX, iy += INCY ) Y[iy] += ALPHA * X[ix];
    else      
      for( i = 0; i < N; ++i, ix += INCX, iy += INCY ) Y[iy] += X[ix] * ALPHA;

  }
};




}; // namespace HAXX
