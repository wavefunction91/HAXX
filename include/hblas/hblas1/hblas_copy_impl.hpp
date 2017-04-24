/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_COPY_IMPL_HPP
#define __INCLUDED_HBLAS_COPY_IMPL_HPP

#include "hblas/hblas1_def.hpp"

namespace HAXX {

/**
 *  Copies the elememts from one quaternion vector to another
 *
 *  Written by DBWY (4/2017)
 *
 *  Based on the BLAS implementaion of DCOPY by Jack Dongarra
 *    http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
 */
template <typename _F>
void HBLAS_COPY(HAXX_INT N, quaternion<_F> *X, HAXX_INT INCX, 
  quaternion<_F> *Y, HAXX_INT INCY) {

  if( N <= 0 ) return;
  // FIXME: See further comments on negative stride
  assert(INCX > 0);
  assert(INCY > 0);

  HAXX_INT i;

  if( INCX == 1 and INCY == 1 ) {

    HAXX_INT m = N % HAXX_COPY_UNROLL;
    if( m != 0 ) {
      for( i = 0; i < m; i++ ) Y[i] = X[i];
      if( N < HAXX_COPY_UNROLL ) return;
    }

    // FIXME: This assumes HAXX_COPY_UNROLL = 7
    for( i = m; i < N; i += HAXX_COPY_UNROLL ) {
      Y[i]   = X[i];
      Y[i+1] = X[i+1];
      Y[i+2] = X[i+2];
      Y[i+3] = X[i+3];
      Y[i+4] = X[i+4];
      Y[i+5] = X[i+5];
      Y[i+6] = X[i+6];
    }

  } else {

    HAXX_INT ix(0), iy(0);
    // FIXME: the original _COPY function has code here to handle
    //   negative increments. Unsure on what that accomplishes

    for( i = 0; i < N; ++i, ix += INCX, iy += INCY ) Y[iy] = X[ix];
      

  }
};

}; // namespace HAXX


#endif
