/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HLAPACK_LARFG_IMPL_HPP
#define __INCLUDED_HLAPACK_LARFG_IMPL_HPP

#include "hlapack/householder_def.hpp"

namespace HAXX {

/**
 *  Written by DBWY (5/2017)
 *
 *  Naming convention and interface based on LAPACK implementation
 *  of ZLARGF. Implementation based in part on Algorithm A1 of
 *
 *    Bunse-Gerstner, A.; Byers, R.; Mehrmann, V.; 
 *      ``A Quaternion QR Algorithm"; Numer. Math. 55, 83-95 (1989)
 */
template < typename _F >
_F HLAPACK_LARFG(HAXX_INT N, quaternion<_F> *X, HAXX_INT INCX){

  _F mu(0.);
  for( HAXX_INT h = 0; h < N; ++h ) mu += norm(X[h]);

  _F ALPHA(0.);
  for( HAXX_INT h = 0; h < N; ++h ) {
    X[h] = X[h] / mu;
    ALPHA = ALPHA + norm(X[h])*norm(X[h]);
  }
  ALPHA = std::sqrt(ALPHA);
  _F TAU = 1./(ALPHA * (ALPHA + norm(X[0])));

  quaternion<_F> sigma;
  // FIXME: Need a full check on X[0] == 0 here
  if( X[0].real() != 0. ) sigma = X[0] / norm(X[0]);
  else             sigma = 1.;

  X[0] += sigma * ALPHA;


/*
quaternion<_F> HLAPACK_LARFG(HAXX_INT N, quaternion<_F>& ALPHA, 
  quaternion<_F> *X, HAXX_INT INCX){

  _F XNORM = HBLAS_NRM2(N-1,X,INCX); // ||x2||
  _F ANORM = norm(ALPHA);            // |ALPHA|

  quaternion<_F> BETA = 
    -std::copysign(std::sqrt(XNORM*XNORM + ANORM*ANORM),ALPHA.real());
  quaternion<_F> TAU = (BETA.real() - ALPHA / BETA.real());
  ALPHA = inv(ALPHA - BETA.real());
  HBLAS_SCAL('L',N-1,ALPHA,X,INCX);
  ALPHA = BETA;
*/

  return TAU;
};

}; // namspace HAXX
#endif
