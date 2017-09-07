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


  // FORTRAN HBLAS1 functions
  extern "C" {
  
    void hscald_(const char*, const HAXX_INT*, const double*,
      const quaternion<double> *, const HAXX_INT*);
    void hscalc_(const char*, const HAXX_INT*, const std::complex<double> *,
      const quaternion<double> *, const HAXX_INT*);
    void hscalh_(const char*, const HAXX_INT*, const quaternion<double> *,
      const quaternion<double> *, const HAXX_INT*);
  
    void hdotu_(const quaternion<double> *, const HAXX_INT*, 
      const quaternion<double> *, const HAXX_INT*, const quaternion<double> *, 
      const HAXX_INT*);
    void hdotc_(const quaternion<double> *, const HAXX_INT*, 
      const quaternion<double> *, const HAXX_INT*, const quaternion<double> *, 
      const HAXX_INT*);
  
    void haxpydh_(const char*, const HAXX_INT*, const double*,
      const quaternion<double> *, const HAXX_INT*, const quaternion<double> *, 
      const HAXX_INT*);
    void haxpych_(const char*, const HAXX_INT*, const std::complex<double> *,
      const quaternion<double> *, const HAXX_INT*, const quaternion<double> *, 
      const HAXX_INT*);
    void haxpyhh_(const char*, const HAXX_INT*, const quaternion<double> *,
      const quaternion<double> *, const HAXX_INT*, const quaternion<double>*, 
      const HAXX_INT*);
  };


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
  void HBLAS_SWAP(const HAXX_INT N, quaternion<_F> * const X, 
    const HAXX_INT INCX, quaternion<_F> * const Y, const HAXX_INT INCY);

  /// Scale a quaternion array by a scalar
  template <typename _F, typename _AlphaF>
  void HBLAS_SCAL(const char SIDE, const HAXX_INT N, const _AlphaF ALPHA, 
    quaternion<_F> * const X, const HAXX_INT INCX);

  /// Copy a quaternion array to another quaternion array
  template <typename _F>
  void HBLAS_COPY(const HAXX_INT N, quaternion<_F> * const X, 
    const HAXX_INT INCX, quaternion<_F> * const Y, const HAXX_INT INCY);

  /// Scale a quaternion array and add it to another quaternion array
  template <typename _F, typename _XF, typename _AlphaF>
  void HBLAS_AXPY(const char SIDE, const HAXX_INT N, const _AlphaF ALPHA, 
    _XF * const X, const HAXX_INT INCX, quaternion<_F> * const Y, 
    const HAXX_INT INCY);

  /// Perform an unaltered dot product of two quaternion arrays
  template <typename _F>
  quaternion<_F> HBLAS_DOTU(const HAXX_INT N, quaternion<_F> * const X, 
    const HAXX_INT INCX, quaternion<_F> * const Y, const HAXX_INT INCY);

  /// Obtain the inner product of two quaternion arrays
  template <typename _F>
  quaternion<_F> HBLAS_DOTC(const HAXX_INT N, quaternion<_F> * const X, 
    const HAXX_INT INCX, quaternion<_F> * const Y, const HAXX_INT INCY);

  /* @} */ // HBLAS1

  /* @} */ // HBLAS

}; // namespace HAXX

#endif
