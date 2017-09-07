/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_DOTC_IMPL_HPP
#define __INCLUDED_HBLAS_DOTC_IMPL_HPP

#include "hblas/hblas1_def.hpp"

namespace HAXX {

/**
 *  Written by DBWY (4/2017)
 *
 *  Based on the BLAS implementation of ZDOTU by Jack Dongarra
 *    http://www.netlib.org/lapack/explore-html/d6/db8/zdotc_8f.html
 *
 * \f$ r,x,y \in \mathbb{H}, \qquad r = \sum_i x^*_i y_i \f$
 */
template <typename _F>
quaternion<_F> HBLAS_DOTC(const HAXX_INT N, quaternion<_F> * const X, 
  const HAXX_INT INCX, quaternion<_F> * const Y, const HAXX_INT INCY) {

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

template<>
quaternion<double> HBLAS_DOTC( HAXX_INT N, quaternion<double> *X, 
  HAXX_INT INCX, quaternion<double> *Y, HAXX_INT INCY){
 
  quaternion<double> htemp;
  hdotc_(&htemp, &N, X, &INCX, Y, &INCY);
  return htemp;

}

}; // namespace HAXX

#endif
