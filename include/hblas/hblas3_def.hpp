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
#include "hblas/hblas_config.hpp"


#define GEMM_FORTRAN_DECL(NAME,F,AMATF,BMATF,ALPHAF,BETAF) \
void NAME##_(const char*, const char*, const HAXX_INT*, const HAXX_INT*,\
  const HAXX_INT*, const ALPHAF*, const AMATF*, const HAXX_INT*,\
  const BMATF*, const HAXX_INT*, const BETAF*, const quaternion<F>*,\
  const HAXX_INT*);


namespace HAXX {

  // FORTRAN HBLAS3 Functions
  extern "C" {
  

    // GEMM functions
      
    GEMM_FORTRAN_DECL(hgemmdd,double,quaternion<double>,quaternion<double>,
      double,double);
    GEMM_FORTRAN_DECL(hgemmdz,double,quaternion<double>,quaternion<double>,
      double,std::complex<double>);
    GEMM_FORTRAN_DECL(hgemmdh,double,quaternion<double>,quaternion<double>,
      double,quaternion<double>);

    GEMM_FORTRAN_DECL(hgemmzd,double,quaternion<double>,quaternion<double>,
      std::complex<double>,double);
    GEMM_FORTRAN_DECL(hgemmzz,double,quaternion<double>,quaternion<double>,
      std::complex<double>,std::complex<double>);
    GEMM_FORTRAN_DECL(hgemmzh,double,quaternion<double>,quaternion<double>,
      std::complex<double>,quaternion<double>);

    GEMM_FORTRAN_DECL(hgemmhd,double,quaternion<double>,quaternion<double>,
      quaternion<double>,double);
    GEMM_FORTRAN_DECL(hgemmhz,double,quaternion<double>,quaternion<double>,
      quaternion<double>,std::complex<double>);
    GEMM_FORTRAN_DECL(hgemmhh,double,quaternion<double>,quaternion<double>,
      quaternion<double>,quaternion<double>);

  };

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
  template <typename _F, typename _AMatF, typename _BMatF, typename _AlphaF, 
    typename _BetaF>
  void HBLAS_GEMM(const char TRANSA, const char TRANSB, const HAXX_INT M, 
    const HAXX_INT N, const HAXX_INT K, const _AlphaF ALPHA, _AMatF * const A, 
    const HAXX_INT LDA, _BMatF * const B, const HAXX_INT LDB, 
    const _BetaF BETA, quaternion<_F> * const C, const HAXX_INT LDC);

  /* @} */ // HBLAS3

  /* @} */ // HBLAS

}; // namespace HAXX

#endif
