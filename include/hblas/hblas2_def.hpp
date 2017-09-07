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
#include "hblas/hblas_config.hpp"

#define GEMV_FORTRAN_DECL(NAME,F,MATF,VECF,ALPHAF,BETAF) \
void NAME##_(const char*, const HAXX_INT*, const HAXX_INT*, const ALPHAF*,\
  const MATF*, const HAXX_INT*, const VECF*, const HAXX_INT*, const BETAF*, \
  const quaternion<F>*, const HAXX_INT*);

#define GER_FORTRAN_DECL(NAME,F,LEFTF,RIGHTF,ALPHAF) \
void NAME##_(const HAXX_INT*, const HAXX_INT*, const ALPHAF*,\
  const LEFTF*, const HAXX_INT*, const RIGHTF*, const HAXX_INT*, \
  const quaternion<F>*, const HAXX_INT*);\

namespace HAXX {

  // FORTRAN HBLAS2 functions
  extern "C" {
  
    // GEMV functions

    GEMV_FORTRAN_DECL(hgemvdd,double,quaternion<double>,quaternion<double>,
      double, double);
    GEMV_FORTRAN_DECL(hgemvdz,double,quaternion<double>,quaternion<double>,
      double, std::complex<double>);
    GEMV_FORTRAN_DECL(hgemvdh,double,quaternion<double>,quaternion<double>,
      double, quaternion<double>);

    GEMV_FORTRAN_DECL(hgemvzd,double,quaternion<double>,quaternion<double>,
      std::complex<double>,double);
    GEMV_FORTRAN_DECL(hgemvzz,double,quaternion<double>,quaternion<double>,
      std::complex<double>,std::complex<double>);
    GEMV_FORTRAN_DECL(hgemvzh,double,quaternion<double>,quaternion<double>,
      std::complex<double>, quaternion<double>);

    GEMV_FORTRAN_DECL(hgemvhd,double,quaternion<double>,quaternion<double>,
      quaternion<double>,double);
    GEMV_FORTRAN_DECL(hgemvhz,double,quaternion<double>,quaternion<double>,
      quaternion<double>,std::complex<double>);
    GEMV_FORTRAN_DECL(hgemvhh,double,quaternion<double>,quaternion<double>,
      quaternion<double>, quaternion<double>);


 
    // GERU functions 
      
    GER_FORTRAN_DECL(hgerud,double,quaternion<double>,quaternion<double>,
      double)
    GER_FORTRAN_DECL(hgeruz,double,quaternion<double>,quaternion<double>,
      std::complex<double>)
    GER_FORTRAN_DECL(hgeruh,double,quaternion<double>,quaternion<double>,
      quaternion<double>)

    // GERC functions

    GER_FORTRAN_DECL(hgercd,double,quaternion<double>,quaternion<double>,
      double)
    GER_FORTRAN_DECL(hgercz,double,quaternion<double>,quaternion<double>,
      std::complex<double>)
    GER_FORTRAN_DECL(hgerch,double,quaternion<double>,quaternion<double>,
      quaternion<double>)

  };

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
  void HBLAS_GEMV(const char TRANS, const HAXX_INT M, const HAXX_INT N, 
    const _AlphaF ALPHA, _MatF * const A, const HAXX_INT LDA, _VecF * const X, 
    const HAXX_INT INCX, const _BetaF BETA, quaternion<_F> * const Y, 
    const HAXX_INT INCY);

  /// Perform the quaternion rank 1 operation with two general vectors
  template <typename _F, typename _LeftVecF, typename _RightVecF, 
    typename _AlphaF>
  void HBLAS_GERU(const HAXX_INT M, const HAXX_INT N, const _AlphaF ALPHA, 
    _LeftVecF * const X, const HAXX_INT INCX, _RightVecF * const Y, 
    const HAXX_INT INCY, quaternion<_F> * const A, const HAXX_INT LDA);

  /// Perform the quaternion rank 1 operation with two general vectors
  template <typename _F, typename _LeftVecF, typename _RightVecF, 
    typename _AlphaF>
  void HBLAS_GERC(const HAXX_INT M, const HAXX_INT N, const _AlphaF ALPHA, 
    _LeftVecF * const X, const HAXX_INT INCX, _RightVecF * const Y, 
    const HAXX_INT INCY, quaternion<_F> * const A, const HAXX_INT LDA);


  /* @} */ // HBLAS2

  /* @} */ // HBLAS
}; // namespace HAXX

#endif
