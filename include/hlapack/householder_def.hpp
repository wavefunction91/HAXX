/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HLAPACK_HOUSEHOLDER_HPP
#define __INCLUDED_HLAPACK_HOUSEHOLDER_HPP

#include "haxx/haxx_def.hpp"
#include "hblas.hpp"

namespace HAXX {

  /**
   *  \addtogroup HLAPACK
   *  @{
   *
   *
   *  @defgroup HLAPACK_HOUSEHOLDER HLAPACK Householder Routines
   *  Householder reflections and manipulations and LAPACK-like algorithms
   *  which utilize them over the quatnerion numbers
   *
   *  @{
   */


  /// Generates a quaternion elementary (Householder) reflector H of order N
  template < typename _F >
/*
  quaternion<_F> HLAPACK_LARFG(HAXX_INT N, quaternion<_F>& ALPHA, 
    quaternion<_F> *X, HAXX_INT INCX);
*/
  _F HLAPACK_LARFG(HAXX_INT N, quaternion<_F> *X, HAXX_INT INCX);

  /// Applies a quaternion elementary (Householder) reflector H of order N
  /// to a matrix C
  template < typename _F >
  void HLAPACK_LARF(char SIDE, HAXX_INT M, HAXX_INT N, quaternion<_F> *V,
    HAXX_INT INCV, quaternion<_F> TAU, quaternion<_F> *C, HAXX_INT LDC,
    quaternion<_F> *WORK);

  template < typename _F >
  HAXX_INT HLAPACK_GEQR2(HAXX_INT M, HAXX_INT N, quaternion<_F> *A, 
    HAXX_INT LDA, quaternion<_F> *TAU, quaternion<_F> *WORK);

  template < typename _F >
  HAXX_INT HLAPACK_GEHD2(HAXX_INT N, HAXX_INT ILO, HAXX_INT IHI,
    quaternion<_F> *A, HAXX_INT LDA, quaternion<_F> *TAU,
    quaternion<_F> *WORK);

  /* @} */ // HLAPACK_HOUSEHOLDER

  /* @} */ // HLAPACK
}; // namespace HAXX
#endif