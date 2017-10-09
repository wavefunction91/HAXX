/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_UTIL_EXPAND_IMPL_HPP
#define __INCLUDED_HBLAS_UTIL_EXPAND_IMPL_HPP

#include <cassert>
#include "hblas/hblas_util.hpp"

namespace HAXX {

  template <>
  void HBLAS_COMPLEX_EXPAND(char ORDER, HAXX_INT M, HAXX_INT N, 
    quaternion<double> *A, HAXX_INT LDA, std::complex<double> *B, 
    HAXX_INT LDB) {

    assert( ORDER == 'F' or ORDER == 'S' );

    if( ORDER == 'F' )
      hzexp1_(&M,&N,A,&LDA,B,&LDB);
    else
      hzexp2_(&M,&N,A,&LDA,B,&LDB);

  }

  template <>
  void HBLAS_REAL_EXPAND(HAXX_INT M, HAXX_INT N, quaternion<double> *A,
    HAXX_INT LDA, double *B, HAXX_INT LDB) {

    hdexp_(&M,&N,A,&LDA,B,&LDB);

  }

}; // namespace HAXX

#endif
