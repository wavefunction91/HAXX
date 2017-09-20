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
 *  Swaps the elements of two strided quaternion arrays of length N
 *
 *  Written by DBWY (4/2017)
 *
 *  Based on the BLAS implementaion of DSWAP by Jack Dongarra
 *    http://www.netlib.org/lapack/explore-html/db/dd4/dswap_8f.html
 */
template <typename _F>
void HBLAS_SWAPV(const HAXX_INT N, quaternion<_F> * const X, 
  const HAXX_INT INCX, quaternion<_F> * const Y, const HAXX_INT INCY) {

  if( N <= 0 ) return;
  // FIXME: See further comments on negative stride
  assert(INCX > 0);
  assert(INCY > 0);

  HAXX_INT i;
  quaternion<_F> qtmp;

  if( INCX == 1 and INCY == 1 ) {
    HAXX_INT m = N % HAXX_SWAP_UNROLL;

    // Use unrolled loops for both unit increments

    // XXX: For some reason Z/CSWAP does not do this. 
    //   Cache utilization?
    if( m != 0 ) {
      for( i = 0; i < m; ++i ) {
        qtmp = X[i];
        X[i] = Y[i];
        Y[i] = qtmp;
      }
      if( N < HAXX_SWAP_UNROLL ) return;
    }

    // FIXME: This assumes HAXX_SWAP_UNROLL = 3
    for( i = m; i < N; i += HAXX_SWAP_UNROLL ) {
      qtmp = X[i];
      X[i] = Y[i];
      Y[i] = qtmp;

      qtmp   = X[i+1];
      X[i+1] = Y[i+1];
      Y[i+1] = qtmp;

      qtmp   = X[i+2];
      X[i+2] = Y[i+2];
      Y[i+2] = qtmp;
    }

  } else {

    HAXX_INT ix(0), iy(0);
    // FIXME: the original _SWAP function has code here to handle
    //   negative increments. Unsure on what that accomplishes

    for( i = 0; i < N; ++i, ix += INCX, iy += INCY ) { 
      qtmp = X[ix];
      X[ix] = Y[iy];
      Y[iy] = qtmp;
    }
  }

};

}; // namespace HAXX


