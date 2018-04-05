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
#include "util/types.hpp"
#include "util/simd.hpp"

namespace HAXX {


  /**
   *  \brief Generic implementation of packing operations used in
   *  the packing utilities.
   *
   *  Will work for an arbitrary data type
   */ 
  template <typename T, typename _TypeWrapper = GenericType<T> >
  struct GenericPackOps {
  
    typedef _TypeWrapper TypeWrapper;
  
    typedef typename _TypeWrapper::load_t     load_t;
    typedef typename _TypeWrapper::noscalar_t noscalar_t;
  

    /**
     *  Generic implementation of scaling operations used prior to internal
     *  operations.
     */
    template <typename U>
    static load_t preOP(load_t &x, U &alpha){ return Mul(alpha,x); } 


    template <typename U, typename... Args>
    static load_t Load(U &alpha, Args... args) {

      auto x = TypeWrapper::Load(args...);
      return preOP(x,alpha);

    }
  
    /**
     *  Generic implementation of scaling operations used prior to internal
     *  operations. Case when no scaling is needed.
     */
  
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
  template < typename T, typename _TypeWrapper = GenericType<T>,
    template<typename,typename> class PackOps = GenericPackOps >
  struct ConjPackOps : public PackOps<T,_TypeWrapper> {
  
    typedef typename  PackOps<T,_TypeWrapper>::load_t load_t;
  
    template <typename U>
    static load_t preOP(load_t &x, U &alpha){ 
      auto y = _TypeWrapper::Conj(x); 
      return PackOps<T,_TypeWrapper>::preOP(y,alpha); 
    }
  
    template <typename U, typename... Args>
    static load_t Load(U &alpha, Args... args) {

      auto x = _TypeWrapper::Load(args...);
      return preOP(x,alpha);

    }

  };
  
  
  
  /**
   *  \brief Specialization of GenericPackOps for double precision
   *  quaternions on AVX / AVX2
   */ 
  template<>
  struct GenericPackOps< quaternion<double>, AVXType<double> > {
  
    typedef quaternion<double>    qd;
    typedef AVXType<double>       TypeWrapper;
  
    typedef typename TypeWrapper::noscalar_t noscalar_t;
    typedef typename TypeWrapper::real_t real_t;
    typedef typename TypeWrapper::complex_t complex_t;
    typedef typename TypeWrapper::quaternion_t quaternion_t;
    typedef typename TypeWrapper::load_t     load_t;
  
  
    template <typename U>
    static inline __m256d preOP(__m256d &x, U &z){ return TypeWrapper::Mul(z,x); }
  
    template <typename U, typename... Args>
    static load_t Load(U &alpha, Args... args) {

      auto x = TypeWrapper::Load(args...);
      return preOP(x,alpha);

    }
  
  
  };
  


  // Forward decl of specialized packing operations

  
  template < typename T = quaternion<double>, 
    typename _TypeWrapper = AVXType<double> > struct GenericPackOps_T1;
  
  template < typename T = quaternion<double>, 
    typename _TypeWrapper = AVXType<double> > struct GenericPackOps_T2;
  
  
  template <typename T = quaternion<double>, 
    typename _TypeWrapper = AVXType<double>>
  using ConjPackOps_T1 = ConjPackOps<T,_TypeWrapper,GenericPackOps_T1>;
  
  template <typename T = quaternion<double>, 
    typename _TypeWrapper = AVXType<double>>
  using ConjPackOps_T2 = ConjPackOps<T,_TypeWrapper,GenericPackOps_T2>;
  
  
  
  
  /**
   *  Factor expensive SIMD permutation into packing for LHS
   *  of quaternion--quaternion matrix product.
   *
   *  Specializes the GenericPackOps for quaternion<double> on
   *  AVX / AVX2
   */ 
  template<>
  struct GenericPackOps_T1< quaternion<double>, AVXType<double> >: 
    public GenericPackOps< quaternion<double>, AVXType<double> > {
  
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
  struct GenericPackOps_T2< quaternion<double>, AVXType<double> >: 
    public GenericPackOps< quaternion<double>, AVXType<double> > {
  
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
