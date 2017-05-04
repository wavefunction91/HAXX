/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HLAPACK_GEHD2_IMPL_HPP
#define __INCLUDED_HLAPACK_GEHD2_IMPL_HPP

#include "hlapack/householder_def.hpp"

namespace HAXX {

/**
 *  Written by DBWY (5/2017)
 *
 *  Naming convention and interface based on LAPACK implementation
 *  of ZGEHD2. Implementation based in part on Algorithm A3 of
 *
 *    Bunse-Gerstner, A.; Byers, R.; Mehrmann, V.; 
 *      ``A Quaternion QR Algorithm"; Numer. Math. 55, 83-95 (1989)
 *
 *  This routine departs from the original Bunse-Gerstner implementation
 *  in that it make use of HLAPACK_LARF application of the reflector
 *  from the right as well as the left.
 */
template < typename _F >
HAXX_INT HLAPACK_GEHD2(HAXX_INT N, HAXX_INT ILO, HAXX_INT IHI,
  quaternion<_F> *A, HAXX_INT LDA, quaternion<_F> *TAU,
  quaternion<_F> *WORK){

  for( HAXX_INT I = ILO; I < IHI - 1; ++I ) {
    TAU[I] = HLAPACK_LARFG(IHI-I,A + std::min(I+2,N) + I*LDA,1);
    HLAPACK_LARF('R',IHI,IHI-I,A + I+1 + I*LDA,1,TAU[I],
      A + (I+1)*LDA,LDA,WORK);
    HLAPACK_LARF('L',IHI-I,N-I,A + I+1 + I*LDA,1,TAU[I],
      A + I+1 (I+1)*LDA,LDA,WORK);
  }

}

}; // namespace HAXX

#endif
