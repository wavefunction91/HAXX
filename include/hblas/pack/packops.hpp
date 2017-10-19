/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_PACK_PACKOPS_HPP
#define __INCLUDED_HBLAS_PACK_PACKOPS_HPP

#include "haxx.hpp"
#include "util/simd.hpp"

namespace HAXX {


  /**
   *  \brief Generic implementation of packing operations used in
   *  the packing utilities.
   *
   *  Will work for an arbitrary data type
   */ 
  template <typename T, typename _TypeWrapper = GenericTypeWrapper<T> >
  struct GenericPackOps {
  
    typedef _TypeWrapper TypeWrapper;
  
    typedef typename _TypeWrapper::load_t     load_t;
    typedef typename _TypeWrapper::noscalar_t noscalar_t;
  

    /// Generic implementation of memory load
    static load_t load(T* x){ return load_t(*x); }

    /**
     * Accounts for fringe case where the dimension being packed
     * is not divisible by the packing size. Loads zero
     */
    static load_t load()    { return load_t(0.); }
  
    /// Generic implementation of store
    static void store(T* x, load_t y){ *x = y; }
  
  
    /**
     *  Generic implementation of scaling operations used prior to internal
     *  operations.
     */
    template <typename U>
    static load_t preOP(load_t &x, U &alpha){ return alpha * x; } 
  
    /**
     *  Generic implementation of scaling operations used prior to internal
     *  operations. Case when no scaling is needed.
     */
    static load_t preOP(load_t &x, noscalar_t &z) { return x; }
  
    /// Generic implementation of (no-op) internal packing operation.
    template <class Tuple>
    static Tuple OP(Tuple &t){ return t; }
  
  };
  
  
  /**
   *  \brief Generic wrapper around a packing implementation which factors
   *  a conjugation operation into the preOP function prior to scaling.
   *
   *  The passed TypeWrapper must have conjOp defined.
   */ 
  template < typename T, typename _TypeWrapper = GenericTypeWrapper<T>,
    template<typename,typename> class PackOps = GenericPackOps >
  struct ConjPackOps : public PackOps<T,_TypeWrapper> {
  
    typedef typename  PackOps<T,_TypeWrapper>::load_t load_t;
  
    template <typename U>
    static load_t preOP(load_t &x, U &alpha){ 
      auto y = _TypeWrapper::conjOp(x); 
      return PackOps<T,_TypeWrapper>::preOP(y,alpha); 
    }
  
  };
  
  
  
  /**
   *  \brief Specialization of GenericPackOps for double precision
   *  quaternions on AVX / AVX2
   */ 
  template<>
  struct GenericPackOps< quaternion<double>, AVX64BitTypeWrapper > {
  
    typedef quaternion<double>    qd;
    typedef AVX64BitTypeWrapper   TypeWrapper;
  
    typedef typename TypeWrapper::noscalar_t noscalar_t;
    typedef typename TypeWrapper::real_t real_t;
    typedef typename TypeWrapper::complex_t complex_t;
    typedef typename TypeWrapper::quaternion_t quaternion_t;
    typedef typename TypeWrapper::load_t     load_t;
  
  
    /// Load quaternion as __m256d
    static inline __m256d load(qd *x){ 
      return LOAD_256D_UNALIGNED_AS(double,x);
    }
  
    /// Fringe case (see GenericPackOps<T>), bcast 0 to vector
    static inline __m256d load(){ return _mm256_setzero_pd(); }
  
    /// Store quaternion as __m256d
    static inline void store(qd* x, __m256d &y) {
      STORE_256D_ALIGNED_AS(double,x,y);
    }
  
  
    /// (No-op) scaling operation
    static inline __m256d preOP(__m256d &x, noscalar_t &z){ return x; }
  
    /// Real scaling operation
    static inline __m256d preOP(__m256d &x, real_t     &z){
      return _mm256_mul_pd(z.x,x);
    }
  
    /// Complex scaling operation
    static inline __m256d preOP(__m256d &x, complex_t     &z){
      __m256d p1 = _mm256_mul_pd(x,z.x);
      __m256d p2 = _mm256_mul_pd(x,z.y);
  
      return _mm256_hsub_pd(p1,p2);
    }
  
  
  
  };
  


  // Forward decl of specialized packing operations

  
  template < typename T = quaternion<double>, 
    typename _TypeWrapper = AVX64BitTypeWrapper > struct GenericPackOps_T1;
  
  template < typename T = quaternion<double>, 
    typename _TypeWrapper = AVX64BitTypeWrapper > struct GenericPackOps_T2;
  
  
  template <typename T = quaternion<double>, 
    typename _TypeWrapper = AVX64BitTypeWrapper>
  using ConjPackOps_T1 = ConjPackOps<T,_TypeWrapper,GenericPackOps_T1>;
  
  template <typename T = quaternion<double>, 
    typename _TypeWrapper = AVX64BitTypeWrapper>
  using ConjPackOps_T2 = ConjPackOps<T,_TypeWrapper,GenericPackOps_T2>;
  
  
  
  
  /**
   *  Factor expensive SIMD permutation into packing for LHS
   *  of quaternion--quaternion matrix product.
   *
   *  Specializes the GenericPackOps for quaternion<double> on
   *  AVX / AVX2
   */ 
  template<>
  struct GenericPackOps_T1< quaternion<double>, AVX64BitTypeWrapper >: 
    public GenericPackOps< quaternion<double>, AVX64BitTypeWrapper > {
  
    using twoTuple = std::tuple<__m256d,__m256d>;
  
    static inline twoTuple  OP( twoTuple &t ) {
      return std::make_tuple(
        _mm256_permute2f128_pd(std::get<0>(t),std::get<1>(t), 0x20),
        _mm256_permute2f128_pd(std::get<0>(t),std::get<1>(t), 0x31)
      );
    }
  
  };
  
  
  
  
  /**
   *  Factor expensive SIMD unpacking into packing for RHS
   *  of quaternion--quaternion matrix product.
   *
   *  Specializes the GenericPackOps for quaternion<double> on
   *  AVX / AVX2
   */ 
  template<>
  struct GenericPackOps_T2< quaternion<double>, AVX64BitTypeWrapper >: 
    public GenericPackOps< quaternion<double>, AVX64BitTypeWrapper > {
  
    using twoTuple = std::tuple<__m256d,__m256d>;
  
    static inline twoTuple  OP( twoTuple &t ) {
      return std::make_tuple(
        _mm256_unpacklo_pd(std::get<0>(t),std::get<1>(t)),
        _mm256_unpackhi_pd(std::get<0>(t),std::get<1>(t))
      );
    }
  
  };


}; // namespace HAXX

#endif
