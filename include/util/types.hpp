/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_UTIL_TYPES_HPP__
#define __INCLUDED_UTIL_TYPES_HPP__

namespace HAXX {


  struct null_t {}; // Null struct to handle various situations


  template <typename T, typename _LoadT = T>
  struct GenericType {

    typedef T      value_type;  // value_type
    typedef _LoadT load_t;      // Load typename
    typedef null_t noscalar_t;  // Null load type

    // Ratio of load_t / T
    static constexpr size_t sz_ratio = sizeof(load_t) / sizeof(T);

    // Sanity checks
    static_assert( 
      sizeof(load_t) % sizeof(T) == 0, 
      "sizeof(load_t) must be an integer multiple of sizeof(T)"
    );

    static_assert( 
      sizeof(load_t) >= sizeof(T), 
      "sizeof(load_t) must >= sizeof(T)"
    );


    // Load from memory
    static inline load_t Load( const T* const x ){ return load_t(*x); }
    static inline load_t Load(            ){ return load_t(0.); }

    // Store to memory
    static inline void Store( T* const x, const load_t y ){ *x = y; }

    // Conjugate
    static inline load_t     Conj( const load_t     &x ){ return SmartConj(x); }
    static inline noscalar_t Conj( const noscalar_t &x ){ return noscalar_t{}; }

    // Load a passed scalar without modification
    template <typename U>
    static inline U LoadScalar( const U &alpha ){ return alpha; }

    // LoadScalar with no scalar Passed
    static inline noscalar_t LoadScalar(){ return noscalar_t{}; }






    // Multiply load_t with anything
    template <typename U, typename V>
    static inline load_t Mul(const U &a, const V &b){ return a * b; }

    // Specializations for noscalar_t
    static inline load_t Mul( const load_t &a, const noscalar_t &b ){ 

      return a;

    }

    static inline load_t Mul( const noscalar_t &a, const load_t &b ){

      return b;

    }



  }; // GenericType



  // Forward decl

  template <typename T> struct AVXType;

  template<>
  struct AVXType<double> : public GenericType<double, __m256d> {

    // Scalar types
    struct real_t       { __m256d x;            }; ///< Real (double)
    struct complex_t    { __m256d x; __m256d y; }; ///< Complex (double)
    struct quaternion_t { __m256d x;            }; ///< Quaternion (double)
    
    // Aliases
    using d  = double;
    using cd = std::complex<double>;
    using qd = quaternion<double>;





    // Inherit
    using GenericType<double,__m256d>::LoadScalar;
    using GenericType<double,__m256d>::Mul;

    // Load from memory
    static inline __m256d Load( const qd * const x ){ 

      return LOAD_256D_UNALIGNED_AS(d,x); 

    }
    static inline __m256d Load(){ return _mm256_setzero_pd(); }


    // Store to memory
    static inline void Store( qd * const x, const __m256d &y ){

      STORE_256D_UNALIGNED_AS(d,x,y);

    }



    // Conjugate
    static inline __m256d Conj( const __m256d &x){ return QCONJ_256D(x); }


    // Load double as q = { x, x, x, x }
    static inline real_t LoadScalar( const double &x ){

      return real_t{ _mm256_broadcast_sd(&x) };

    }


    // Load complex double as q = { { x(R), x(I),x(R), x(I) }, 
    //                              { x(I),-x(R),x(I),-x(R) } }
    static inline complex_t LoadScalar( const std::complex<double> & x ){

      const __m256i maskConj = _mm256_set_epi64x(
                                 0x8000000000000000, 0,
                                 0x8000000000000000, 0 );
  
      __m128d alphaC = LOAD_128D_UNALIGNED_AS(double,&x);
      __m256d alpha  = SET_256D_FROM_128D(alphaC,alphaC);
      __m256d alpha_C =
        _mm256_xor_pd(_mm256_permute_pd(alpha, 0x5),
                      _mm256_castsi256_pd(maskConj)
                     );
  
      return complex_t{ alpha, alpha_C };

    }

    // Load quaternion double as is
    static inline quaternion_t LoadScalar( const qd &x ){

      return quaternion_t{ Load(&x) };

    }




    // Multiply load_t with real_t
    static inline __m256d Mul( const __m256d &a, const real_t &b ) {

      return _mm256_mul_pd(b.x, a);

    }
    static inline __m256d Mul( const real_t &a, const __m256d &b ) {

      return _mm256_mul_pd(a.x, b);

    }

    // Multiply load_t with complex_t
    static inline __m256d Mul( const complex_t &a, const __m256d &b ) {

      __m256d p1 = _mm256_mul_pd(b, a.x);
      __m256d p2 = _mm256_mul_pd(b, a.y);

      return _mm256_hsub_pd(p1,p2);

    }
    static inline __m256d Mul( const __m256d &a, const complex_t &b ) {

      // FIXME: NEEDS VALIDATAION
      const __m256i mask1 = _mm256_set_epi64x(0x8000000000000000,0,0,0);
      const __m256i mask2 = _mm256_set_epi64x(0,0x8000000000000000,0,0);


      __m256d b1 = _mm256_xor_pd(b.x,_mm256_castsi256_pd(mask1));
      __m256d b2 = _mm256_xor_pd(b.y,_mm256_castsi256_pd(mask2));

      __m256d p1 = _mm256_mul_pd(a, b1);
      __m256d p2 = _mm256_mul_pd(a, b2);

      return _mm256_hsub_pd(p1,p2);

    }

    // Multiply load_t with quaternion_t
    static inline __m256d Mul( const __m256d &a, const quaternion_t &b ){

      return MULDQ_NN(a,b.x);

    }
    static inline __m256d Mul( const quaternion_t &a, const __m256d &b ){

      return MULDQ_NN(a.x,b);

    }


  }; // AVXType<double>


};


#endif
