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

  /// Multiply a general vector by a quaternion matrix
  template <typename _F,typename _MatF, typename _VecF, typename _AlphaF, 
    typename _BetaF>
  void GEMV(char TRANS, HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, 
    _MatF *A, HAXX_INT LDA, _VecF *X, HAXX_INT INCX, 
    _BetaF BETA, quaternion<_F> *Y, HAXX_INT INCY);

  /// Perform the quaternion rank 1 operation with two general vectors
  template <typename _F, typename _LeftVecF, typename _RightVecF, typename _AlphaF>
  void GERU(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, _LeftVecF *X,
    HAXX_INT INCX, _RightVecF *Y, HAXX_INT INCY, quaternion<_F> *A, 
    HAXX_INT LDA);

  /// Perform the quaternion rank 1 operation with two general vectors
  template <typename _F, typename _LeftVecF, typename _RightVecF, typename _AlphaF>
  void GERC(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, _LeftVecF *X,
    HAXX_INT INCX, _RightVecF *Y, HAXX_INT INCY, quaternion<_F> *A, 
    HAXX_INT LDA);


  /* @} */ // HBLAS2

  /* @} */ // HBLAS
}; // namespace HAXX

#endif
