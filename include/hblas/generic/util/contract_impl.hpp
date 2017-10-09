/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_UTIL_CONTRACT_IMPL_HPP
#define __INCLUDED_HBLAS_UTIL_CONTRACT_IMPL_HPP

#include "hblas/hblas_util.hpp"

namespace HAXX {

  template <>
  void HBLAS_COMPLEX_CONTRACT(char ORDER, char UPLO, HAXX_INT M, HAXX_INT N, 
    quaternion<double> *A, HAXX_INT LDA, std::complex<double> *B, 
    HAXX_INT LDB) {

    assert( ORDER == 'F' or ORDER == 'S' );

    if( ORDER == 'F' )
      hzcon1_(&UPLO,&M,&N,A,&LDA,B,&LDB);
    else
      hzcon2_(&UPLO,&M,&N,A,&LDA,B,&LDB);

  }


}; // namespace HAXX

#endif
