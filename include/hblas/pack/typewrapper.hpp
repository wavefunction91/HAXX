/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_PACK_TYPEWRAPPER_HPP
#define __INCLUDED_HBLAS_PACK_TYPEWRAPPER_HPP

#include "haxx.hpp"
#include "util/simd.hpp"

namespace HAXX {


  /**
   *  \brief Generic implementation of scalar operations used in
   *  the packing utilities. 
   *
   *  Will work with any scalar data type.
   */
  template <typename T>
  struct GenericTypeWrapper {
  
    typedef T load_t; ///< Load variable typename
  
    /// Struct to handle the case when no scalar is passed to
    /// the packing utility
    struct noscalar_t {};
  
  
    /// Load the passed scalar without modification
    template <typename U>
    static const U cacheScalar( U &alpha ){ return alpha; }
  
    /// No scalar is passed
    static noscalar_t cacheScalar(){ return noscalar_t{}; }
  

    /// Generic, typesafe implementation of conjugation
    template <typename U>
    static U conjOp( U &t ) { return SmartConj(t); }
  
  };
  

  /**
   *  \brief AVX / AVX2 (256-bit vector length) implementation
   *  of scalar operations used in packing utilities.
   *
   *  Currently only viable for quaternion packing (FIXME)
   */ 
  struct AVX64BitTypeWrapper {
  
    typedef __m256d load_t; ///< Load variable typename
  
    struct noscalar_t   {};                        ///< No scalar passed
    struct real_t       { __m256d x;            }; ///< Real (double)
    struct complex_t    { __m256d x; __m256d y; }; ///< Complex (double)
    struct quaternion_t { __m256d x;            }; ///< Quaternion (double)
  
  
  
    /// No scalar passed to packing utility
    static inline noscalar_t cacheScalar(){ return noscalar_t(); }
  
    /// Real (double) scalar passed to packing utility (bcast to __m256d)
    static inline const real_t cacheScalar( double &alpha ) {
      return real_t{_mm256_broadcast_sd(&alpha)};
    }
  
    /**
     *  Complex (double) scalar passed to packing utility. 
     *  Load as
     *    {
     *      { x, x },
     *      { -x*, -x* }
     *    }
     *  
     */ 
    static inline const complex_t 
      cacheScalar(std::complex<double> &ALPHA) {
  
      const __m256i maskConj = _mm256_set_epi64x(
                                 0x8000000000000000, 0,
                                 0x8000000000000000, 0 );
  
      __m128d alphaC = LOAD_128D_UNALIGNED_AS(double,&ALPHA);
      __m256d alpha  = SET_256D_FROM_128D(alphaC,alphaC);
      __m256d alpha_C =
        _mm256_xor_pd(_mm256_permute_pd(alpha, 0x5),
                      _mm256_castsi256_pd(maskConj)
                     );
  
      return complex_t{ alpha, alpha_C };
    }
  
    /// Quaternion (double) scalar passed to packing utility (load as __m256d)
    static inline const quaternion_t cacheScalar( quaternion<double>& alpha ) {
      return quaternion_t{LOAD_256D_UNALIGNED_AS(double,&alpha)};
    }

    /// AVX / AVX2 Conjugation operation (FIXME: only works for quaternions) 
    static __m256d conjOp( __m256d &t ) { return QCONJ_256D(t); }
  
  };

};

#endif
