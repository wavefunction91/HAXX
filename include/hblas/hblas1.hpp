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

// Preprocessor macros for quick FORTRAN declarations

#define SCAL_FORTRAN_DECL(NAME,F,ALPHAF) \
void NAME##_(const char*, const HAXX_INT*, const ALPHAF*, \
  const quaternion<F> *, const HAXX_INT*);

#define DOT_FORTRAN_DECL(NAME,F) \
void NAME##_(const quaternion<F> *, const HAXX_INT*, const quaternion<F> *,\
  const HAXX_INT*, const quaternion<F> *, const HAXX_INT*);

#define AXPY_FORTRAN_DECL(NAME,F,XF,ALPHAF)\
void NAME##_(const char*, const HAXX_INT*, const ALPHAF*, const XF *, \
  const HAXX_INT*, const quaternion<F> *, const HAXX_INT*);





namespace HAXX {


  // FORTRAN HBLAS1 functions
  extern "C" {

    // SCAL functions  
    SCAL_FORTRAN_DECL(hscald,double,double);
    SCAL_FORTRAN_DECL(hscalc,double,std::complex<double>);
    SCAL_FORTRAN_DECL(hscalh,double,quaternion<double>);

    // DOT functions
    DOT_FORTRAN_DECL(hdotu,double);
    DOT_FORTRAN_DECL(hdotc,double);

    // AXPY functions
    AXPY_FORTRAN_DECL(haxpydh,double,quaternion<double>,double);
    AXPY_FORTRAN_DECL(haxpych,double,quaternion<double>,std::complex<double>);
    AXPY_FORTRAN_DECL(haxpyhh,double,quaternion<double>,quaternion<double>);

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


  /**
   *  @defgroup HBLAS1V Level 1 HBLASV
   *  Level 1 BLAS operations on quaternion vectors
   *
   *  @{
   */ 

  /// Swap the states of two quaternion arrays
  template <typename _F>
  void HBLAS_SWAPV(const HAXX_INT N, quaternion<_F> * const X, 
    const HAXX_INT INCX, quaternion<_F> * const Y, const HAXX_INT INCY);

  /// Copy a quaternion array to another quaternion array
  template <typename _F>
  void HBLAS_COPYV(const HAXX_INT N, quaternion<_F> * const X, 
    const HAXX_INT INCX, quaternion<_F> * const Y, const HAXX_INT INCY);

  /// Scale a quaternion array by a scalar
  template <typename _F, typename _AlphaF>
  void HBLAS_SCALV(const char SIDE, const HAXX_INT N, const _AlphaF ALPHA, 
    quaternion<_F> * const X, const HAXX_INT INCX);

  /// Scale a quaternion array and add it to another quaternion array
  template <typename _F, typename _XF, typename _AlphaF>
  void HBLAS_AXPYV(const char SIDE, const HAXX_INT N, const _AlphaF ALPHA, 
    _XF * const X, const HAXX_INT INCX, quaternion<_F> * const Y, 
    const HAXX_INT INCY);

  /// Perform an unaltered dot product of two quaternion arrays
  template <typename _F>
  quaternion<_F> HBLAS_DOTUV(const HAXX_INT N, quaternion<_F> * const X, 
    const HAXX_INT INCX, quaternion<_F> * const Y, const HAXX_INT INCY);

  /// Obtain the inner product of two quaternion arrays
  template <typename _F>
  quaternion<_F> HBLAS_DOTCV(const HAXX_INT N, quaternion<_F> * const X, 
    const HAXX_INT INCX, quaternion<_F> * const Y, const HAXX_INT INCY);

  /* @} */ // HBLAS1V


  /**
   *  @defgroup HBLAS1M Level 1 HBLASM
   *  Level 1 BLAS operations on quaternion matricies (ala BLIS)
   *
   *  @{
   */ 

  /// Scale / Transpose / Conjugate a quaternion matrix in place
  template <typename _F, typename _AlphaF>
  void HBLAS_SCALM(const char SIDE, const char TRANSA, const HAXX_INT M, 
    const HAXX_INT N, const _AlphaF ALPHA, quaternion<_F> * const A, 
    const HAXX_INT LDA, const HAXX_INT INCA);

  /// Accumulate a quaternion matrix
  template <typename _F, typename _XF, typename _AlphaF>
  void HBLAS_AXPYM(const char SIDE, const char TRANSA, const HAXX_INT M, 
    const HAXX_INT N, const _AlphaF ALPHA, _XF * const A, const HAXX_INT LDA, 
    const HAXX_INT INCA, quaternion<_F> * const B, const HAXX_INT LDB, 
    const HAXX_INT INCB);

  /* @} */ // HBLAS1M

  /* @} */ // HBLAS1

  /* @} */ // HBLAS

}; // namespace HAXX

#endif
