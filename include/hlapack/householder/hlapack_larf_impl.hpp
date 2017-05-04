/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HLAPACK_LARF_IMPL_HPP
#define __INCLUDED_HLAPACK_LARF_IMPL_HPP

#include "hlapack/householder_def.hpp"

namespace HAXX {

/**
 *  Written by DBWY (5/2017)
 *
 *  Naming convention and interface based on LAPACK implementation
 *  of ZLARF. Implementation based in part on Algorithm A2 of
 *
 *    Bunse-Gerstner, A.; Byers, R.; Mehrmann, V.; 
 *      ``A Quaternion QR Algorithm"; Numer. Math. 55, 83-95 (1989)
 */
template < typename _F >
void HLAPACK_LARF(char SIDE, HAXX_INT M, HAXX_INT N, quaternion<_F> *V,
  HAXX_INT INCV, quaternion<_F> TAU, quaternion<_F> *C, HAXX_INT LDC,
  quaternion<_F> *WORK) {

  if( SIDE == 'L' ) {
    // FIXME: This is stupid, should generalize GEMV to multiply the
    // vector from the left and to allow for conjugation of X. 
    //   This implementaion assumes that INCV == 1
    assert( INCV == 1 );

    HBLAS_GEMM('C','N',1,N,M,_F(1.),V,M,C,LDC,_F(0.),WORK,1);
    HBLAS_GERU(M,N,-TAU,V,INCV,WORK,1,C,LDC);
  } else {

    HBLAS_GEMV('N',M,N,_F(1.),C,LDC,V,INCV,_F(0.),WORK,1);
    // FIXME: This assumes that TAU is real
    HBLAS_GERC(M,N,-TAU,WORK,1,V,INCV,C,LDC);

  }
}

}; // namespace HAXX
#endif
