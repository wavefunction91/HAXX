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
quaternion<_F> HLAPACK_LARFG(HAXX_INT N, quaternion<_F>& ALPHA, 
  quaternion<_F> *X, HAXX_INT INCX){

/*
  double mu(0.);
  for( HAXX_INT h = 0; h < N; ++h ) mu += norm(X[h]);

  quaternion<_F> TAU(0.);
  for( HAXX_INT h = 0; h < N; ++h ) {
    X[h] = X[h] / mu;
    ALPHA = ALPHA + norm(X[h])*norm(X[h]);
  }
  ALPHA = sqrt(ALPHA);
  TAU = inv(ALPHA * (ALPHA + norm(X[h])));

  quaternion<_F> sigma;
  if( X[0] != 0. ) sigma = X[0] / norm(X[0]);
  else             sigma = 1.;

  X[0] += sigma * ALPHA;

  return TAU;
*/



};

}; // namspace HAXX
#endif
