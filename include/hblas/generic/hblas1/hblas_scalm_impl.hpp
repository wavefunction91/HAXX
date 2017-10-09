/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#pragma once

#include "hblas/hblas1.hpp"

namespace HAXX {

/**
 *  Scale / Transpose / Conjugate a quaternion matrix in place
 *
 *  Written by DBWY (9/2017)
 */
template <typename _F, typename _AlphaF>
void HBLAS_SCALM(const char SIDE, const char TRANSA, const HAXX_INT M, 
  const HAXX_INT N, const _AlphaF ALPHA, quaternion<_F> * const A, 
  const HAXX_INT LDA, const HAXX_INT INCA) {

  if( N <= 0 or M <= 0 or INCA <= 0 ) return;
  assert(TRANSA == 'N'); // Only supporting scaling for now

  HAXX_INT j;

  quaternion<_F> *locA = A;

  for( j = 0; j < N; j++ ) {
    HBLAS_SCALV(SIDE,M,ALPHA,locA,INCA);
    locA += LDA;
  }


};



}; // namespace HAXX

