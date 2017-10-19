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
#include "util/boilerplate.hpp"



namespace HAXX {







template <typename T>
struct GenericTypeWrapper {

  typedef T load_t;

  struct noscalar_t {};


  template <typename U>
  static const U cacheScalar( U &alpha ){ return alpha; }

  static noscalar_t cacheScalar(){ return noscalar_t{}; }

  template <typename U>
  static U conjOp( U &t ) { return SmartConj(t); }

};

struct AVX64BitTypeWrapper {

  typedef __m256d load_t;

  struct noscalar_t   {};
  struct real_t       { __m256d x;            };
  struct complex_t    { __m256d x; __m256d y; };
  struct quaternion_t { __m256d x;            };



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

  static inline const quaternion_t cacheScalar( quaternion<double>& alpha ) {
    return quaternion_t{LOAD_256D_UNALIGNED_AS(double,&alpha)};
  }

  static __m256d conjOp( __m256d &t ) { return QCONJ_256D(t); }

};




template <typename T, typename _TypeWrapper = GenericTypeWrapper<T> >
struct GenericPackOps {

  typedef _TypeWrapper TypeWrapper;

  typedef typename _TypeWrapper::load_t     load_t;
  typedef typename _TypeWrapper::noscalar_t noscalar_t;

  static load_t load(T* x){ return load_t(*x); }
  static load_t load()    { return load_t(0.); }

  static void store(T* x, load_t y){ *x = y; }


  template <typename U>
  static load_t preOP(load_t &x, U &alpha){ return alpha * x; } 

  static load_t preOP(load_t &x, noscalar_t &z) { return x; }

  template <class Tuple>
  static Tuple OP(Tuple &t){ return t; }

};


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



template<>
struct GenericPackOps< quaternion<double>, AVX64BitTypeWrapper > {

  typedef quaternion<double>    qd;
  typedef AVX64BitTypeWrapper   TypeWrapper;

  typedef typename TypeWrapper::noscalar_t noscalar_t;
  typedef typename TypeWrapper::real_t real_t;
  typedef typename TypeWrapper::complex_t complex_t;
  typedef typename TypeWrapper::quaternion_t quaternion_t;
  typedef typename TypeWrapper::load_t     load_t;


  static inline __m256d load(qd *x){ 
    return LOAD_256D_UNALIGNED_AS(double,x);
  }

  static inline __m256d load(){ return _mm256_setzero_pd(); }

  static inline void store(qd* x, __m256d &y) {
    STORE_256D_ALIGNED_AS(double,x,y);
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





template <size_t PK, typename T, typename PackOps, typename... Args>
void TPACK(const HAXX_INT M, const HAXX_INT N, T *X,
  const HAXX_INT LDX, T *Xp, Args... args) {


  // Sanity Check
  static_assert(
    (sizeof(typename PackOps::load_t) % sizeof(T)) == 0,
    "The size of the loaded variable must be an integer multiple of T" 
  );



  HAXX_INT i,j;


  T *_x,*_x1,*_x2;
  T *_xp;
  

  _x  = X;
  _xp = Xp;

  auto alpha = PackOps::TypeWrapper::cacheScalar(args...);


  auto load_preOP = 
    [&](auto &t){ 
      auto x = PackOps::load(t);        
      return PackOps::preOP(x,alpha);
    };

  auto ptrInc = [](auto &t){ return t + 1; };

  auto store = [](auto *&p, auto &t){ PackOps::store(p,t);    };

  j = N / PK;
  if( j > 0 )
  do {

    auto xCol = ptrSeq<PK>(_x,LDX);

    _x += PK*LDX;

    for( i = 0; i < M; i++ ) {

      auto x   = apply( load_preOP, xCol );
      auto x_t = PackOps::OP(x);

      arrExc( store, _xp, x_t );

      xCol = apply( ptrInc, xCol );
      _xp += PK;

    }
    
    j--;

  } while(j > 0);

  if( N % PK ) {

    #ifdef _NDEBUG
      #undef _NDEBUG
      assert(false /* This has yet to be worked out */);
      #define _NDEBUG
    #endif

#if 0
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
#endif

  }

}

// NPACK
template <size_t PK, typename T, typename PackOps, typename... Args>
inline void NPACK(const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp, Args... args) {

  // Sanity Check
  static_assert(
    (sizeof(typename PackOps::load_t) % sizeof(T)) == 0,
    "The size of the loaded variable must be an integer multiple of T" 
  );

  HAXX_INT i,j;
  T *_x,*_x1,*_x2;
  T *_xp;

  _x = X;
  _xp = Xp;

  auto alpha = PackOps::TypeWrapper::cacheScalar(args...);

  auto load_preOP = 
    [&](auto &t){ 
      auto x = PackOps::load(t);        
      return PackOps::preOP(x,alpha);
    };

  auto ptrInc = [&](auto &t){ return t + LDX; };

  auto store = [](auto *&p, auto &t){ PackOps::store(p,t);    };

  i = M / PK;
  if( i > 0 )
  do {

    auto xCol = ptrSeq<PK>(_x,1);

    _x += PK;

    for( j = 0; j < N; j++ ) {

      auto x   = apply( load_preOP, xCol );
      auto x_t = PackOps::OP(x);
      arrExc( store, _xp, x_t );

      xCol = apply( ptrInc, xCol );
      _xp += PK;

    }

    i--;
  } while( i > 0 );

  if( M % PK ) {

    #ifdef _NDEBUG
      #undef _NDEBUG
      assert(false /* This has yet to be worked out */);
      #define _NDEBUG
    #endif

#if 0
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
#endif

  }

};



}; // namespace HAXX
