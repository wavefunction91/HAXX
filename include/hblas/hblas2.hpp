/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS2_HPP
#define __INCLUDED_HBLAS2_HPP

#include "haxx/haxx_def.hpp"

#ifndef HAXX_INT
  #define HAXX_INT int
#endif

#define RANK2_INDX(i,j,N) ( (i) + (j)*(N) )

namespace HAXX {

  /**
   *  \addtogroup HBLAS
   *  @{
   *
   *
   *  @defgroup HBLAS2 Level 2 HBLAS
   *  Level 2 BLAS operations over quaternion numbers
   *
   *  @{
   */

  /// Multiply a quaternion vector by a quaternion matrix
  template <typename _F, typename _AlphaF, typename _BetaF>
  void GEMV(char TRANS, HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, 
    quaternion<_F> *A, HAXX_INT LDA, quaternion<_F> *X, 
    HAXX_INT INCX, _BetaF BETA, quaternion<_F> *Y, HAXX_INT INCY);

  /// Multiply a complex vector by a quaternion matrix
  template <typename _F, typename _AlphaF, typename _BetaF>
  void GEMV(char TRANS, HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, 
    quaternion<_F> *A, HAXX_INT LDA, std::complex<_F> *X, 
    HAXX_INT INCX, _BetaF BETA, quaternion<_F> *Y, HAXX_INT INCY);

  /// Multiply a real vector by a quaternion matrix
  template <typename _F, typename _AlphaF, typename _BetaF>
  void GEMV(char TRANS, HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, 
    quaternion<_F> *A, HAXX_INT LDA, _F *X, HAXX_INT INCX,
    _BetaF BETA, quaternion<_F> *Y, HAXX_INT INCY);



  /// Perform the rank 1 operation with quaternion vectors
  template <typename _F, typename _AlphaF>
  void GERU(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, quaternion<_F> *X,
    HAXX_INT INCX, quaternion<_F> *Y, HAXX_INT INCY, quaternion<_F> *A, 
    HAXX_INT LDA);

  /// Perform the rank 1 operation with one quaternion and one complex vector
  template <typename _F, typename _AlphaF>
  void GERU(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, quaternion<_F> *X,
    HAXX_INT INCX, std::complex<_F> *Y, HAXX_INT INCY, quaternion<_F> *A, 
    HAXX_INT LDA);
  /// Perform the rank 1 operation with one quaternion and one complex vector
  template <typename _F, typename _AlphaF>
  void GERU(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, std::complex<_F> *X, 
    HAXX_INT INCX, quaternion<_F> *Y, HAXX_INT INCY, quaternion<_F> *A, 
    HAXX_INT LDA);

  /// Perform the rank 1 operation with one quaternion and one real vector
  template <typename _F, typename _AlphaF>
  void GERU(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, quaternion<_F> *X,
    HAXX_INT INCX, _F *Y, HAXX_INT INCY, quaternion<_F> *A, HAXX_INT LDA);
  /// Perform the rank 1 operation with one quaternion and one real vector
  template <typename _F, typename _AlphaF>
  void GERU(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, _F *X, HAXX_INT INCX, 
    quaternion<_F> *Y, HAXX_INT INCY, quaternion<_F> *A, HAXX_INT LDA);

  /// Perform the rank 1 operation with quaternion vectors
  template <typename _F, typename _AlphaF>
  void GERC(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, quaternion<_F> *X,
    HAXX_INT INCX, quaternion<_F> *Y, HAXX_INT INCY, quaternion<_F> *A, 
    HAXX_INT LDA);

  /// Perform the rank 1 operation with one quaternion and one complex vector
  template <typename _F, typename _AlphaF>
  void GERC(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, quaternion<_F> *X,
    HAXX_INT INCX, std::complex<_F> *Y, HAXX_INT INCY, quaternion<_F> *A, 
    HAXX_INT LDA);
  /// Perform the rank 1 operation with one quaternion and one complex vector
  template <typename _F, typename _AlphaF>
  void GERC(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, std::complex<_F> *X, 
    HAXX_INT INCX, quaternion<_F> *Y, HAXX_INT INCY, quaternion<_F> *A, 
    HAXX_INT LDA);

  /// Perform the rank 1 operation with one quaternion and one real vector
  template <typename _F, typename _AlphaF>
  void GERC(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, quaternion<_F> *X,
    HAXX_INT INCX, _F *Y, HAXX_INT INCY, quaternion<_F> *A, HAXX_INT LDA);
  /// Perform the rank 1 operation with one quaternion and one real vector
  template <typename _F, typename _AlphaF>
  void GERC(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, _F *X, HAXX_INT INCX, 
    quaternion<_F> *Y, HAXX_INT INCY, quaternion<_F> *A, HAXX_INT LDA);

  /* @} */ // HBLAS2

  /* @} */ // HBLAS
}; // namespace HAXX

#endif
