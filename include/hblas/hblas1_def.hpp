/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS1_HPP
#define __INCLUDED_HBLAS1_HPP

#include "haxx/haxx_def.hpp"
#include "hblas/hblas_config.hpp"

namespace HAXX {

  /**
   *  \addtogroup HBLAS
   *  @{
   *
   *
   *  @defgroup HBLAS1 Level 1 HBLAS
   *  Level 1 BLAS operations over quaternion numbers
   *
   *  @{
   */

  /// Swap the states of two quaternion arrays
  template <typename _F>
  void HBLAS_SWAP(HAXX_INT N, quaternion<_F> *X, HAXX_INT INCX, 
    quaternion<_F> *Y, HAXX_INT INCY);

  /// Scale a quaternion array by a scalar
  template <typename _F, typename _AlphaF>
  void HBLAS_SCAL(char SIDE, HAXX_INT N, _AlphaF ALPHA, quaternion<_F> *X, 
    HAXX_INT INCX);

  /// Copy a quaternion array to another quaternion array
  template <typename _F>
  void HBLAS_COPY(HAXX_INT N, quaternion<_F> *X, HAXX_INT INCX, 
    quaternion<_F> *Y, HAXX_INT INCY);

  /// Scale a quaternion array and add it to another quaternion array
  template <typename _F, typename _XF, typename _AlphaF>
  void HBLAS_AXPY(char SIDE, HAXX_INT N, _AlphaF ALPHA, _XF *X, 
    HAXX_INT INCX, quaternion<_F> *Y, HAXX_INT INCY);

  /// Perform an unaltered dot product of two quaternion arrays
  template <typename _F>
  quaternion<_F> HBLAS_DOTU( HAXX_INT N, quaternion<_F> *X, HAXX_INT INCX, 
    quaternion<_F> *Y, HAXX_INT INCY);

  /// Obtain the inner product of two quaternion arrays
  template <typename _F>
  quaternion<_F> HBLAS_DOTC( HAXX_INT N, quaternion<_F> *X, HAXX_INT INCX, 
    quaternion<_F> *Y, HAXX_INT INCY);

  /* @} */ // HBLAS1

  /* @} */ // HBLAS

}; // namespace HAXX

#endif
