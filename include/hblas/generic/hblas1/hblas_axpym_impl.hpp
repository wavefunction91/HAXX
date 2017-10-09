/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#pragma once

#include <cassert>
#include "hblas/hblas1.hpp"

namespace HAXX {

/**
 *  Written by DBWY (4/2017)
 */
template <typename _F, typename _XF, typename _AlphaF>
void HBLAS_AXPYM(const char SIDE, const char TRANSA, const HAXX_INT M, 
  const HAXX_INT N, const _AlphaF ALPHA, _XF * const A, const HAXX_INT LDA, 
  const HAXX_INT INCA, quaternion<_F> * const B, const HAXX_INT LDB, 
  const HAXX_INT INCB) {

  
  if( N <= 0 or M <= 0) return;
  if( ALPHA == _AlphaF(0.) ) return;

  assert(TRANSA == 'N'); // Only supporting scaling for now

  HAXX_INT j;
  
  quaternion<_F> *locA = A, *locB = B;

  for( j = 0; j < N; j++ ) {
    HBLAS_AXPYV(SIDE,M,ALPHA,locA,INCA,locB,INCB);
    locA += LDA;
    locB += LDB;
  }

};




}; // namespace HAXX
