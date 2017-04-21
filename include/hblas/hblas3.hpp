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

#ifndef HAXX_INT
  #define HAXX_INT int
#endif

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
  template <typename _F, typename _AlphaF, typename _BetaF>
  void GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT, N, HAXX_INT K,
    _AlphaF ALPHA, quaternion<_F> *A, HAXX_INT LDA, quaternion<_F> *B, 
    HAXX_INT LDB, _BetaF BETA, quaternion<_F> *C, HAXX_INT LDC);

  /// Multiply a quaternion matrix by a complex matrix
  template <typename _F, typename _AlphaF, typename _BetaF>
  void GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT, N, HAXX_INT K,
    _AlphaF ALPHA, quaternion<_F> *A, HAXX_INT LDA, std::complex<_F> *B, 
    HAXX_INT LDB, _BetaF BETA, quaternion<_F> *C, HAXX_INT LDC);

  /// Multiply a complex matrix by a quaternion matrix
  template <typename _F, typename _AlphaF, typename _BetaF>
  void GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT, N, HAXX_INT K,
    _AlphaF ALPHA, std::complex<_F> *A, HAXX_INT LDA, quaternion<_F> *B, 
    HAXX_INT LDB, _BetaF BETA, quaternion<_F> *C, HAXX_INT LDC);

  /// Multiply a quaternion matrix by a real matrix
  template <typename _F, typename _AlphaF, typename _BetaF>
  void GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT, N, HAXX_INT K,
    _AlphaF ALPHA, quaternion<_F> *A, HAXX_INT LDA, _F *B, 
    HAXX_INT LDB, _BetaF BETA, quaternion<_F> *C, HAXX_INT LDC);

  /// Multiply a real matrix by a quaternion matrix
  template <typename _F, typename _AlphaF, typename _BetaF>
  void GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT, N, HAXX_INT K,
    _AlphaF ALPHA, _F *A, HAXX_INT LDA, quaternion<_F> *B, 
    HAXX_INT LDB, _BetaF BETA, quaternion<_F> *C, HAXX_INT LDC);

  /* @} */ // HBLAS3

  /* @} */ // HBLAS

}; // namespace HAXX
