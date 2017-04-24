/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_SCAL_IMPL_HPP
#define __INCLUDED_HBLAS_SCAL_IMPL_HPP

#include "hblas/hblas1_def.hpp"

namespace HAXX {

/**
 *  Scales a quaternion vector in place
 *
 *  Written by DBWY (4/2017)
 *
 *  Based on the BLAS implementaion of DSCAL by Jack Dongarra
 *    http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html 
 *
 *  SIDE == 'L'
 *
 *  \f$ x_i = \alpha x_i \f$
 *
 *  SIDE == 'R'
 *
 *  \f$ x_i = x_i \alpha  \f$
 */
template <typename _F, typename _AlphaF>
void HBLAS_SCAL(char SIDE, HAXX_INT N, _AlphaF ALPHA, quaternion<_F> *X, 
  HAXX_INT INCX){

  if( N <= 0 or INCX <= 0 ) return;

  HAXX_INT i;

  // FIXME: Should write a specialization for real ALPHA where side
  //   doesnt matter
  bool isR = SIDE == 'R';
  bool isL = not isR;

  if( INCX == 1 ) {
    HAXX_INT m = N % HAXX_SCAL_UNROLL;
  
    if( m != 0 ) {
      if( isL ) for( i = 0; i < m; ++i ) X[i] = ALPHA * X[i];
      else      for( i = 0; i < m; ++i ) X[i] = X[i] * ALPHA;

      if( N < HAXX_SCAL_UNROLL ) return;
    }

    // FIXME: This assumes HAXX_SCAL_UNROLL = 5
    if( isL )
      for( i = m; i < N; i += HAXX_SCAL_UNROLL ) {
        X[i]   = ALPHA * X[i];
        X[i+1] = ALPHA * X[i+1];
        X[i+2] = ALPHA * X[i+2];
        X[i+3] = ALPHA * X[i+3];
        X[i+4] = ALPHA * X[i+4];
      }
    else
      for( i = m; i < N; i += HAXX_SCAL_UNROLL ) {
        X[i]   = X[i]   * ALPHA; 
        X[i+1] = X[i+1] * ALPHA; 
        X[i+2] = X[i+2] * ALPHA; 
        X[i+3] = X[i+3] * ALPHA; 
        X[i+4] = X[i+4] * ALPHA; 
      }
  } else {

    HAXX_INT NINCX = N*INCX;
    if( isL ) for( i = 0; i < NINCX; i += INCX ) X[i] = ALPHA * X[i];
    else      for( i = 0; i < NINCX; i += INCX ) X[i] = X[i] * ALPHA;

  }


};

}; // namespace HAXX

#endif
