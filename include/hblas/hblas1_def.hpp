/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS1_IMPL_HPP
#define __INCLUDED_HBLAS1_IMPL_HPP

#include "hblas/hblas1.hpp"
#include <cassert>
#include <iostream>

// Hardcoded unrolling parameters
#define HAXX_SWAP_UNROLL 3
#define HAXX_SCAL_UNROLL 5
#define HAXX_COPY_UNROLL 7
#define HAXX_AXPY_UNROLL 4

namespace HAXX {


/**
 *  Swaps the elements of two strided quaternion arrays of length N
 *
 *  Based on the BLAS implementaion of DSWAP by Jack Dongarra
 *    http://www.netlib.org/lapack/explore-html/db/dd4/dswap_8f.html
 *
 *
 */
template <typename _F>
void SWAP(HAXX_INT N, quaternion<_F> *X, HAXX_INT INCX, quaternion<_F> *Y, 
  HAXX_INT INCY) {

  if( N <= 0 ) return;
  // FIXME: See further comments on negative stride
  assert(INCX > 0);
  assert(INCY > 0);

  HAXX_INT i;
  quaternion<_F> qtmp;

  if( INCX == 1 and INCY == 1 ) {
    HAXX_INT m = N % HAXX_SWAP_UNROLL;

    // Use unrolled loops for both unit increments

    // FIXME: For some reason Z/CSWAP does not do this. 
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

/**
 *  Scales a quaternion vector in place
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
 *
 *  
 */
template <typename _F, typename _AlphaF>
void SCAL(char SIDE, HAXX_INT N, _AlphaF ALPHA, quaternion<_F> *X, 
  HAXX_INT INCX){

  if( N <= 0 or INCX <= 0 ) return;

  HAXX_INT i;

  // FIXME: Should write a specialization for real ALPHA where side
  //   doesnt matter
  bool isR = SIDE == 'R';
  bool isL = not isR;

  if( INCX == 1 ) {
    // FIXME: This seems strange to hard code this
    //   I would assuming that this is arch dependent for 
    //   best performace
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


/**
 *  Copies the elememts from one quaternion vector to another
 *
 *  Based on the BLAS implementaion of DCOPY by Jack Dongarra
 *    http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
 */
template <typename _F>
void COPY(HAXX_INT N, quaternion<_F> *X, HAXX_INT INCX, 
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
    for( i = 0; i < N; i += HAXX_COPY_UNROLL ) {
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

/**
 *  Based on the BLAS implementaion of DCOPY by Jack Dongarra
 *    http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
 *
 *
 *  SIDE == 'L'
 *
 *  \f$ y_i = \alpha x_i  + y_i\f$
 *
 *  SIDE == 'R'
 *
 *  \f$ y_i = x_i \alpha + y_i \f$
 */
template <typename _F, typename _AlphaF>
void AXPY(char SIDE, HAXX_INT N, _AlphaF ALPHA, quaternion<_F> *X, 
  HAXX_INT INCX, quaternion<_F> *Y, HAXX_INT INCY) {

  
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

      // FIXME: DAXPY has this outside of the if-check? Unline COPY and SCAL
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

/**
 * Based on the BLAS implementation of ZDOTU by Jack Dongarra
 *   http://www.netlib.org/lapack/explore-html/db/d2d/zdotu_8f.html 
 *
 * \f$ r,x,y \in \mathbb{H}, \qquad r = \sum_i x_i y_i
 */
template <typename _F>
quaternion<_F> DOTU( HAXX_INT N, quaternion<_F> *X, HAXX_INT INCX, 
  quaternion<_F> *Y, HAXX_INT INCY){

  quaternion<_F> htemp(0.,0.,0.,0.);
  if( N <= 0 ) return htemp;

  // FIXME: See further comments on negative stride
  assert(INCX > 0);
  assert(INCY > 0);

  HAXX_INT i;

  if( INCX == 1 and INCY == 1 ) {

    for( i = 0; i < N; ++i ) htemp += X[i]*Y[i];

  } else {

    HAXX_INT ix(0), iy(0);
    // FIXME: the original _AXPY function has code here to handle
    //   negative increments. Unsure on what that accomplishes
      
    for( i = 0; i < N; ++i, ix += INCX, iy += INCY ){ 
      htemp += X[ix] * Y[iy];
    }
  }  

  return htemp;
};

/**
 * Based on the BLAS implementation of ZDOTU by Jack Dongarra
 *   http://www.netlib.org/lapack/explore-html/d6/db8/zdotc_8f.html
 *
 * \f$ r,x,y \in \mathbb{H}, \qquad r = \sum_i x^*_i y_i
 */
template <typename _F>
quaternion<_F> DOTC( HAXX_INT N, quaternion<_F> *X, HAXX_INT INCX, 
  quaternion<_F> *Y, HAXX_INT INCY){

  quaternion<_F> htemp(0.,0.,0.,0.);

  if( N <= 0 ) return htemp;

  // FIXME: See further comments on negative stride
  assert(INCX > 0);
  assert(INCY > 0);

  HAXX_INT i;

  if( INCX == 1 and INCY == 1 ) {

    for( i = 0; i < N; ++i ) htemp += conj(X[i])*Y[i];

  } else {

    HAXX_INT ix(0), iy(0);
    // FIXME: the original _AXPY function has code here to handle
    //   negative increments. Unsure on what that accomplishes
      
    for( i = 0; i < N; ++i, ix += INCX, iy += INCY ) 
      htemp += conj(X[ix]) * Y[iy];
  }  

  return htemp;
};

}; // namespace HAXX

#endif
