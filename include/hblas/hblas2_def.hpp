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


namespace HAXX {

  // FORTRAN HBLAS2 functions
  extern "C" {
  
    // GEMV functions
    void hgemvdd_(const char*, const HAXX_INT*, const HAXX_INT*, const double*,
      const quaternion<double>*, const HAXX_INT*, const quaternion<double>*, 
      const HAXX_INT*, const double*, const quaternion<double>*, 
      const HAXX_INT*);

    void hgemvdz_(const char*, const HAXX_INT*, const HAXX_INT*, const double*,
      const quaternion<double>*, const HAXX_INT*, const quaternion<double>*, 
      const HAXX_INT*, const std::complex<double>*, const quaternion<double>*, 
      const HAXX_INT*);

    void hgemvdh_(const char*, const HAXX_INT*, const HAXX_INT*, const double*,
      const quaternion<double>*, const HAXX_INT*, const quaternion<double>*, 
      const HAXX_INT*, const quaternion<double>*, const quaternion<double>*, 
      const HAXX_INT*);

    void hgemvzd_(const char*, const HAXX_INT*, const HAXX_INT*, 
      const std::complex<double>*, const quaternion<double>*, const HAXX_INT*, 
      const quaternion<double>*, const HAXX_INT*, const double*, 
      const quaternion<double>*, const HAXX_INT*);

    void hgemvzz_(const char*, const HAXX_INT*, const HAXX_INT*, 
      const std::complex<double>*, const quaternion<double>*, const HAXX_INT*, 
      const quaternion<double>*, const HAXX_INT*, const std::complex<double>*, 
      const quaternion<double>*, const HAXX_INT*);

    void hgemvzh_(const char*, const HAXX_INT*, const HAXX_INT*, 
      const std::complex<double>*, const quaternion<double>*, const HAXX_INT*, 
      const quaternion<double>*, const HAXX_INT*, const quaternion<double>*, 
      const quaternion<double>*, const HAXX_INT*);

    void hgemvhd_(const char*, const HAXX_INT*, const HAXX_INT*, 
      const quaternion<double>*, const quaternion<double>*, const HAXX_INT*, 
      const quaternion<double>*, const HAXX_INT*, const double*, 
      const quaternion<double>*, const HAXX_INT*);

    void hgemvhz_(const char*, const HAXX_INT*, const HAXX_INT*, 
      const quaternion<double>*, const quaternion<double>*, const HAXX_INT*, 
      const quaternion<double>*, const HAXX_INT*, const std::complex<double>*, 
      const quaternion<double>*, const HAXX_INT*);

    void hgemvhh_(const char*, const HAXX_INT*, const HAXX_INT*, 
      const quaternion<double>*, const quaternion<double>*, const HAXX_INT*, 
      const quaternion<double>*, const HAXX_INT*, const quaternion<double>*, 
      const quaternion<double>*, const HAXX_INT*);


  
  
  
  
    // GERU functions
  
    void hgerud_(const HAXX_INT*, const HAXX_INT*, const double*, 
      const quaternion<double>*, const HAXX_INT*, const quaternion<double>*, 
      const HAXX_INT*, const quaternion<double>*, const HAXX_INT*);

    void hgeruz_(const HAXX_INT*, const HAXX_INT*, const std::complex<double>*, 
      const quaternion<double>*, const HAXX_INT*, const quaternion<double>*, 
      const HAXX_INT*, const quaternion<double>*, const HAXX_INT*);

    void hgeruh_(const HAXX_INT*, const HAXX_INT*, const quaternion<double>*, 
      const quaternion<double>*, const HAXX_INT*, const quaternion<double>*, 
      const HAXX_INT*, const quaternion<double>*, const HAXX_INT*);
  

    // GERC functions
    
    void hgercd_(const HAXX_INT*, const HAXX_INT*, const double*, 
      const quaternion<double>*, const HAXX_INT*, const quaternion<double>*, 
      const HAXX_INT*, const quaternion<double>*, const HAXX_INT*);

    void hgercz_(const HAXX_INT*, const HAXX_INT*, const std::complex<double>*, 
      const quaternion<double>*, const HAXX_INT*, const quaternion<double>*, 
      const HAXX_INT*, const quaternion<double>*, const HAXX_INT*);

    void hgerch_(const HAXX_INT*, const HAXX_INT*, const quaternion<double>*, 
      const quaternion<double>*, const HAXX_INT*, const quaternion<double>*, 
      const HAXX_INT*, const quaternion<double>*, const HAXX_INT*);

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
