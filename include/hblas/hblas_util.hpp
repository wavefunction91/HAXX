/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_UTIL_HPP
#define __INCLUDED_HBLAS_UTIL_HPP

#include "haxx/haxx_def.hpp" 
#include "hblas/hblas_config.hpp"

extern "C" {

  void hzexp1_(const int*, const int*, HAXX::quaternion<double>*, const int*,
    std::complex<double>*, const int*);
  void hzexp2_(const int*, const int*, HAXX::quaternion<double>*, const int*,
    std::complex<double>*, const int*);
  void hdexp_(const int*, const int*, HAXX::quaternion<double>*, const int*,
    double*, const int*);

}

namespace HAXX {

  /**
   *  \addtogroup HBLAS
   *  @{
   *
   *
   *  @defgroup HBLAS_UTIL HBLAS Utilities
   *  Utility functions for HBLAS
   *
   *  @{
   */


  /// Expand a quaternion matrix to a complex matrix
  template <typename _F>
  void HBLAS_COMPLEX_EXPAND(char ORDER, HAXX_INT M, HAXX_INT N, 
    quaternion<_F> *A, HAXX_INT LDA, std::complex<_F> *B, HAXX_INT LDB);

  /// Expans a quaternion matrix to a real matrix
  template <typename _F>
  void HBLAS_REAL_EXPAND(HAXX_INT M, HAXX_INT N, quaternion<_F> *A,
    HAXX_INT LDA, _F *B, HAXX_INT LDB);



  /* @} */ // HBLAS_UTIL

  /* @} */ // HBLAS

}; // namespace HAXX

#endif
