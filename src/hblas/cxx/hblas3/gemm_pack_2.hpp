/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#include "haxx.hpp"
#include "util/simd.hpp"

namespace HAXX {

template <typename T>
struct GenericPackOpsBase {

  struct noscalar_t {};

  static T load(T* x){ return *x; }
  static T load()    { return T(0.); }

  static void store(T* x, T y){ *x = y; }

  template <typename U>
  static const U cacheScalar( U &alpha ){ return alpha; }

  static noscalar_t cacheScalar(){ return noscalar_t(); }



  template <typename U>
  static const T preOP(T &x, U &alpha){ return alpha * x; } 

  static const T preOP(T &x, noscalar_t &z) { return x; }


};

template <typename T>
struct GenericPackOps2 : public GenericPackOpsBase<T> {

  static T OP1(T &x1, T &x2) { return x1; }
  static T OP2(T &x1, T &x2) { return x2; }

};

template <typename T>
struct ConjPackOps2 : public GenericPackOps2<T> {

  template <typename U>
  static const T preOP(T &x, U &alpha){ 
    auto y = SmartConj(x); 
    return GenericPackOpsBase<T>::preOP(y,alpha); 
  } 

};


template<>
struct GenericPackOpsBase< quaternion<double> > {

  typedef quaternion<double> qd;
  struct noscalar_t {};
  struct real_t{ __m256d x; };
  struct complex_t{ __m256d x; __m256d  y; };
  struct quaternion_t{ __m256d x; };

  static inline __m256d load(qd *x){ 
    return LOAD_256D_UNALIGNED_AS(double,x);
  }

  static inline __m256d load(){ return _mm256_setzero_pd(); }

  static inline void store(qd* x, __m256d &y) {
    STORE_256D_ALIGNED_AS(double,x,y);
  }

  static inline noscalar_t cacheScalar(){ return noscalar_t(); }

  static inline const real_t cacheScalar( double &alpha ) {
    return real_t{_mm256_broadcast_sd(&alpha)};
  }

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

  static inline const quaternion_t cacheScalar( qd& alpha ) {
    return quaternion_t{LOAD_256D_UNALIGNED_AS(double,&alpha)};
  }


  static inline __m256d preOP(__m256d &x, noscalar_t &z){ return x; }

  static inline __m256d preOP(__m256d &x, real_t     &z){
    return _mm256_mul_pd(z.x,x);
  }

  static inline __m256d preOP(__m256d &x, complex_t     &z){
    __m256d p1 = _mm256_mul_pd(x,z.x);
    __m256d p2 = _mm256_mul_pd(x,z.y);

    return _mm256_hsub_pd(p1,p2);
  }


};

struct GenericPackOps_T1 : public GenericPackOpsBase< quaternion<double> > {

  static inline __m256d OP1(__m256d &x1, __m256d &x2) {
    return _mm256_permute2f128_pd(x1,x2, 0x20);
  }

  static inline __m256d OP2(__m256d &x1, __m256d &x2) {
    return _mm256_permute2f128_pd(x1,x2, 0x31);
  }

};


struct ConjPackOps_T1 : public GenericPackOps_T1 {

  template <typename U>
  static __m256d preOP(__m256d &x, U &z) {
    auto y = QCONJ_256D(x); return GenericPackOps_T1::preOP(y,z);
  }

};

struct GenericPackOps_T2 : public GenericPackOpsBase< quaternion<double> > {

  static inline __m256d OP1(__m256d &x1, __m256d &x2) {
    return _mm256_unpacklo_pd(x1, x2);
  }

  static inline __m256d OP2(__m256d &x1, __m256d &x2) {
    return _mm256_unpackhi_pd(x1, x2);
  }

};



struct ConjPackOps_T2 : public GenericPackOps_T2 {

  template <typename U>
  static __m256d preOP(__m256d &x, U &z) {
    auto y = QCONJ_256D(x); return GenericPackOps_T2::preOP(y,z);
  }

};




template <typename T, typename PackOps, typename... Args>
void TPACK2(const HAXX_INT M, const HAXX_INT N, T *X,
  const HAXX_INT LDX, T *Xp, Args... args) {


  HAXX_INT i,j;


  T *_x,*_x1,*_x2;
  T *_xp;
  

  _x  = X;
  _xp = Xp;

  auto alpha = PackOps::cacheScalar(args...);

  j = N / 2;

  if( j > 0 )
  do {

    _x1 = _x;
    _x2 = _x1 + LDX;
 
    _x += 2*LDX;

    for( i = 0; i < M; i++ ) {

      auto x1 = PackOps::load(_x1);
      auto x2 = PackOps::load(_x2);

      x1 = PackOps::preOP(x1,alpha);
      x2 = PackOps::preOP(x2,alpha);

      auto x1_t = PackOps::OP1(x1, x2);
      auto x2_t = PackOps::OP2(x1, x2);

      PackOps::store(_xp,   x1_t);
      PackOps::store(_xp+1, x2_t);

      _xp += 2;
      _x1++; _x2++;

    }
    
    j--;

  } while(j > 0);

  if( N % 2 ) {
    for( i = 0; i < M; i++ ) {

      auto x1 = PackOps::load(_x);
      auto x2 = PackOps::load();

      x1 = PackOps::preOP(x1,alpha);
      x2 = PackOps::preOP(x2,alpha);

      auto x1_t = PackOps::OP1(x1, x2);
      auto x2_t = PackOps::OP2(x1, x2);

      PackOps::store(_xp,   x1_t);
      PackOps::store(_xp+1, x2_t);

      _xp +=2; _x++;

    }
  }

}

// NPACK2, no scalar
template <typename T, typename PackOps, typename... Args>
inline void NPACK2(const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp, Args... args) {


  HAXX_INT i,j;
  T *_x,*_x1,*_x2;
  T *_xp;

  _x = X;
  _xp = Xp;

  auto alpha = PackOps::cacheScalar(args...);

  i = M / 2;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 2;


    for( j = 0; j < N; j++ ) {

      auto x1 = PackOps::load(_x1    );
      auto x2 = PackOps::load(_x1 + 1);

      x1 = PackOps::preOP(x1,alpha);
      x2 = PackOps::preOP(x2,alpha);

      auto x1_t = PackOps::OP1(x1, x2);
      auto x2_t = PackOps::OP2(x1, x2);

      PackOps::store(_xp,   x1_t);
      PackOps::store(_xp+1, x2_t);

      _xp += 2; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 2 ) {
    for( j = 0; j < N; j++ ) {

      auto x1 = PackOps::load(_x);
      auto x2 = PackOps::load();

      x1 = PackOps::preOP(x1,alpha);
      x2 = PackOps::preOP(x2,alpha);

      auto x1_t = PackOps::OP1(x1, x2);
      auto x2_t = PackOps::OP2(x1, x2);

      PackOps::store(_xp,   x1_t);
      PackOps::store(_xp+1, x2_t);

      _x += LDX; _xp += 2;

    }
  }

};



}; // namespace HAXX
