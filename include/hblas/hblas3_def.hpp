/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS3_HPP
#define __INCLUDED_HBLAS3_HPP

#include "haxx/haxx_def.hpp"
#include "hblas/hblas_config.hpp"

namespace HAXX {

  /**
   *  \addtogroup HBLAS
   *  @{
   *
   *
   *  @defgroup HBLAS3 Level 3 HBLAS
   *  Level 3 BLAS operations over quaternion numbers
   *
   *  @{
   */

  /// Multiply a quaternion matrix by a quaternion matrix
  template <typename _F, typename _AMatF, typename _BMatF, typename _AlphaF, 
    typename _BetaF>
  void HBLAS_GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT N, HAXX_INT K,
    _AlphaF ALPHA, _AMatF *A, HAXX_INT LDA, _BMatF *B, HAXX_INT LDB, 
    _BetaF BETA, quaternion<_F> *C, HAXX_INT LDC);

  /* @} */ // HBLAS3

  /* @} */ // HBLAS

}; // namespace HAXX

#endif