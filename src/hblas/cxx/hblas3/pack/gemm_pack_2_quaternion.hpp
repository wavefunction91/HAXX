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


// TPACK2, no scalar
template <>
inline void TPACK2(const HAXX_INT M, const HAXX_INT N, 
  quaternion<double> *X, const HAXX_INT LDX, 
  quaternion<double> *Xp) {

  HAXX_INT i,j;


  quaternion<double> *_x,*_x1,*_x2;
  quaternion<double> *_xp;
  

  _x = X;
  _xp = Xp;


  j = N / 2;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
 
    _x += 2*LDX;

    for( i = 0; i < M; i++ ) {

      __m256d x1 = LOAD_256D_UNALIGNED_AS(double,_x1);
      __m256d x2 = LOAD_256D_UNALIGNED_AS(double,_x2);


      STORE_256D_ALIGNED_AS(double,_xp  ,x1);
      STORE_256D_ALIGNED_AS(double,_xp+1,x2);

      _xp += 2;
      _x1++; _x2++;
    }
    
    j--;

  } while(j > 0);

  if( N % 2 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 2); j++) {

        __m256d x = LOAD_256D_UNALIGNED_AS(double,_x+j*LDX);
        STORE_256D_ALIGNED_AS(double,_xp,x);

        _xp++;
      }

      for( j = 0; j < 2 - (N % 2); j++) {
        __m256d x = _mm256_setzero_pd();
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;
      }

      _x++;
    }
  }

};

// TPACK2 with real scalar
template <>
inline void TPACK2(const double ALPHA, const HAXX_INT M, 
  const HAXX_INT N, quaternion<double> *X, 
  const HAXX_INT LDX, quaternion<double> *Xp) {

  HAXX_INT i,j;


  quaternion<double> *_x,*_x1,*_x2;
  quaternion<double> *_xp;
  
  const __m256d alpha = _mm256_broadcast_sd(&ALPHA);

  _x = X;
  _xp = Xp;


  j = N / 2;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
 
    _x += 2*LDX;

    for( i = 0; i < M; i++ ) {

      __m256d x1 = LOAD_256D_UNALIGNED_AS(double,_x1);
      __m256d x2 = LOAD_256D_UNALIGNED_AS(double,_x2);

      x1 = _mm256_mul_pd(alpha,x1);
      x2 = _mm256_mul_pd(alpha,x2);


      STORE_256D_ALIGNED_AS(double,_xp  ,x1);
      STORE_256D_ALIGNED_AS(double,_xp+1,x2);

      _xp += 2;
      _x1++; _x2++;
    }
    
    j--;

  } while(j > 0);

  if( N % 2 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 2); j++) {

        __m256d x = LOAD_256D_UNALIGNED_AS(double,_x+j*LDX);
        x = _mm256_mul_pd(alpha,x);

        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;

      }

      for( j = 0; j < 2 - (N % 2); j++) {
        __m256d x = _mm256_setzero_pd();
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;
      }

      _x++;
    }
  }

};


// TPACK2 with complex scalar
template <>
inline void TPACK2(const std::complex<double> ALPHA, 
  const HAXX_INT M, const HAXX_INT N, quaternion<double> *X, 
  const HAXX_INT LDX, quaternion<double> *Xp) {

  HAXX_INT i,j;


  quaternion<double> *_x,*_x1,*_x2;
  quaternion<double> *_xp;

  const __m256i maskConj = _mm256_set_epi64x(
                             0x8000000000000000, 0,
                             0x8000000000000000, 0 );
  
  const __m128d alphaC = LOAD_128D_UNALIGNED_AS(double,&ALPHA);
  const __m256d alpha  = SET_256D_FROM_128D(alphaC,alphaC);
  const __m256d alpha_conj = 
    _mm256_xor_pd(_mm256_permute_pd(alpha, 0x5),
                  _mm256_castsi256_pd(maskConj)
                 );


  _x = X;
  _xp = Xp;


  j = N / 2;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
 
    _x += 2*LDX;

    for( i = 0; i < M; i++ ) {

      __m256d x1 = LOAD_256D_UNALIGNED_AS(double,_x1);
      __m256d x2 = LOAD_256D_UNALIGNED_AS(double,_x2);

      __m256d p1 = _mm256_mul_pd(x1,alpha);
      __m256d p2 = _mm256_mul_pd(x1,alpha_conj);

      x1 = _mm256_hsub_pd(p1,p2);

      p1 = _mm256_mul_pd(x2,alpha);
      p2 = _mm256_mul_pd(x2,alpha_conj);

      x2 = _mm256_hsub_pd(p1,p2);


      STORE_256D_ALIGNED_AS(double,_xp  ,x1);
      STORE_256D_ALIGNED_AS(double,_xp+1,x2);

      _xp += 2;
      _x1++; _x2++;
    }
    
    j--;

  } while(j > 0);

  if( N % 2 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 2); j++) {

        __m256d x = LOAD_256D_UNALIGNED_AS(double,_x+j*LDX);

        __m256d p1 = _mm256_mul_pd(x,alpha);
        __m256d p2 = _mm256_mul_pd(x,alpha_conj);

        x = _mm256_hsub_pd(p1,p2);

        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;

      }

      for( j = 0; j < 2 - (N % 2); j++) {
        __m256d x = _mm256_setzero_pd();
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;
      }

      _x++;
    }
  }

};


// TPACK2 with conjugation, no scalar
template <>
inline void TPACKC2(const HAXX_INT M, const HAXX_INT N, 
  quaternion<double> *X, const HAXX_INT LDX, 
  quaternion<double> *Xp) {

  HAXX_INT i,j;


  quaternion<double> *_x,*_x1,*_x2;
  quaternion<double> *_xp;
  

  _x = X;
  _xp = Xp;


  j = N / 2;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
 
    _x += 2*LDX;

    for( i = 0; i < M; i++ ) {

      __m256d x1 = LOAD_256D_UNALIGNED_AS(double,_x1);
      __m256d x2 = LOAD_256D_UNALIGNED_AS(double,_x2);

      x1 = QCONJ_256D(x1);
      x2 = QCONJ_256D(x2);

      STORE_256D_ALIGNED_AS(double,_xp  ,x1);
      STORE_256D_ALIGNED_AS(double,_xp+1,x2);

      _xp += 2;
      _x1++; _x2++;
    }
    
    j--;

  } while(j > 0);

  if( N % 2 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 2); j++) {

        __m256d x = LOAD_256D_UNALIGNED_AS(double,_x+j*LDX);

        x = QCONJ_256D(x);

        STORE_256D_ALIGNED_AS(double,_xp,x);

        _xp++;
      }

      for( j = 0; j < 2 - (N % 2); j++) {
        __m256d x = _mm256_setzero_pd();
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;
      }

      _x++;
    }
  }

};


// TPACK2 with conjugation real scalar
template <>
inline void TPACKC2(const double ALPHA, const HAXX_INT M, 
  const HAXX_INT N, quaternion<double> *X, 
  const HAXX_INT LDX, quaternion<double> *Xp) {

  HAXX_INT i,j;


  quaternion<double> *_x,*_x1,*_x2;
  quaternion<double> *_xp;
  
  const __m256d alpha = _mm256_broadcast_sd(&ALPHA);

  _x = X;
  _xp = Xp;


  j = N / 2;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
 
    _x += 2*LDX;

    for( i = 0; i < M; i++ ) {

      __m256d x1 = LOAD_256D_UNALIGNED_AS(double,_x1);
      __m256d x2 = LOAD_256D_UNALIGNED_AS(double,_x2);

      x1 = QCONJ_256D(x1);
      x2 = QCONJ_256D(x2);

      x1 = _mm256_mul_pd(alpha,x1);
      x2 = _mm256_mul_pd(alpha,x2);


      STORE_256D_ALIGNED_AS(double,_xp  ,x1);
      STORE_256D_ALIGNED_AS(double,_xp+1,x2);

      _xp += 2;
      _x1++; _x2++;
    }
    
    j--;

  } while(j > 0);

  if( N % 2 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 2); j++) {

        __m256d x = LOAD_256D_UNALIGNED_AS(double,_x+j*LDX);

        x = QCONJ_256D(x);
        x = _mm256_mul_pd(alpha,x);

        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;

      }

      for( j = 0; j < 2 - (N % 2); j++) {
        __m256d x = _mm256_setzero_pd();
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;
      }

      _x++;
    }
  }

};


// TPACK2 with conjugation and complex scalar
template <>
inline void TPACKC2(const std::complex<double> ALPHA, 
  const HAXX_INT M, const HAXX_INT N, quaternion<double> *X, 
  const HAXX_INT LDX, quaternion<double> *Xp) {

  HAXX_INT i,j;


  quaternion<double> *_x,*_x1,*_x2;
  quaternion<double> *_xp;

  const __m256i maskConj = _mm256_set_epi64x(
                             0x8000000000000000, 0,
                             0x8000000000000000, 0 );
  
  const __m128d alphaC = LOAD_128D_UNALIGNED_AS(double,&ALPHA);
  const __m256d alpha  = SET_256D_FROM_128D(alphaC,alphaC);
  const __m256d alpha_conj = 
    _mm256_xor_pd(_mm256_permute_pd(alpha, 0x5),
                  _mm256_castsi256_pd(maskConj)
                 );


  _x = X;
  _xp = Xp;


  j = N / 2;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
 
    _x += 2*LDX;

    for( i = 0; i < M; i++ ) {

      __m256d x1 = LOAD_256D_UNALIGNED_AS(double,_x1);
      __m256d x2 = LOAD_256D_UNALIGNED_AS(double,_x2);

      x1 = QCONJ_256D(x1);
      x2 = QCONJ_256D(x2);

      __m256d p1 = _mm256_mul_pd(x1,alpha);
      __m256d p2 = _mm256_mul_pd(x1,alpha_conj);

      x1 = _mm256_hsub_pd(p1,p2);

      p1 = _mm256_mul_pd(x2,alpha);
      p2 = _mm256_mul_pd(x2,alpha_conj);

      x2 = _mm256_hsub_pd(p1,p2);


      STORE_256D_ALIGNED_AS(double,_xp  ,x1);
      STORE_256D_ALIGNED_AS(double,_xp+1,x2);

      _xp += 2;
      _x1++; _x2++;
    }
    
    j--;

  } while(j > 0);

  if( N % 2 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 2); j++) {

        __m256d x = LOAD_256D_UNALIGNED_AS(double,_x+j*LDX);

        x = QCONJ_256D(x);

        __m256d p1 = _mm256_mul_pd(x,alpha);
        __m256d p2 = _mm256_mul_pd(x,alpha_conj);

        x = _mm256_hsub_pd(p1,p2);

        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;

      }

      for( j = 0; j < 2 - (N % 2); j++) {
        __m256d x = _mm256_setzero_pd();
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;
      }

      _x++;
    }
  }

};


// NPACK2, no scalar
template <>
inline void NPACK2(const HAXX_INT M, const HAXX_INT N, 
  quaternion<double> *X, const HAXX_INT LDX, 
  quaternion<double> *Xp) {


  HAXX_INT i,j;
  quaternion<double> *_x,*_x1,*_x2;
  quaternion<double> *_xp;

  _x = X;
  _xp = Xp;

  i = M / 2;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 2;


    for( j = 0; j < N; j++ ) {

      __m256d x1 = LOAD_256D_UNALIGNED_AS(double,_x1);
      __m256d x2 = LOAD_256D_UNALIGNED_AS(double,_x1+1);


      STORE_256D_ALIGNED_AS(double,_xp  ,x1);
      STORE_256D_ALIGNED_AS(double,_xp+1,x2);

      _xp += 2; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 2 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 2); i++ ){
        __m256d x = LOAD_256D_UNALIGNED_AS(double,_x+i);
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++; 
      }

      for( i = 0; i < 2 - (M % 2); i++) {
        __m256d x = _mm256_setzero_pd();
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;
      }

      _x += LDX;
    }
  }

};

// NPACK2, with real scalar
template <>
inline void NPACK2(const double ALPHA, const HAXX_INT M, 
  const HAXX_INT N, quaternion<double> *X, 
  const HAXX_INT LDX, quaternion<double> *Xp) {


  HAXX_INT i,j;
  quaternion<double> *_x,*_x1,*_x2;
  quaternion<double> *_xp;

  _x = X;
  _xp = Xp;

  const __m256d alpha = _mm256_broadcast_sd(&ALPHA);

  i = M / 2;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 2;


    for( j = 0; j < N; j++ ) {

      __m256d x1 = LOAD_256D_UNALIGNED_AS(double,_x1);
      __m256d x2 = LOAD_256D_UNALIGNED_AS(double,_x1+1);

      x1 = _mm256_mul_pd(alpha,x1);
      x2 = _mm256_mul_pd(alpha,x2);


      STORE_256D_ALIGNED_AS(double,_xp  ,x1);
      STORE_256D_ALIGNED_AS(double,_xp+1,x2);

      _xp += 2; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 2 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 2); i++ ){
        __m256d x = LOAD_256D_UNALIGNED_AS(double,_x+i);
        x = _mm256_mul_pd(alpha,x);
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++; 
      }

      for( i = 0; i < 2 - (M % 2); i++) {
        __m256d x = _mm256_setzero_pd();
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;
      }

      _x += LDX;
    }
  }

};

// NPACK2, with complex scalar
template <>
inline void NPACK2(const std::complex<double> ALPHA, 
  const HAXX_INT M, const HAXX_INT N, quaternion<double> *X, 
  const HAXX_INT LDX, quaternion<double> *Xp) {


  HAXX_INT i,j;
  quaternion<double> *_x,*_x1,*_x2;
  quaternion<double> *_xp;

  _x = X;
  _xp = Xp;

  const __m256i maskConj = _mm256_set_epi64x(
                             0x8000000000000000, 0,
                             0x8000000000000000, 0 );
  
  const __m128d alphaC = LOAD_128D_UNALIGNED_AS(double,&ALPHA);
  const __m256d alpha  = SET_256D_FROM_128D(alphaC,alphaC);
  const __m256d alpha_conj = 
    _mm256_xor_pd(_mm256_permute_pd(alpha, 0x5),
                  _mm256_castsi256_pd(maskConj)
                 );

  i = M / 2;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 2;


    for( j = 0; j < N; j++ ) {

      __m256d x1 = LOAD_256D_UNALIGNED_AS(double,_x1);
      __m256d x2 = LOAD_256D_UNALIGNED_AS(double,_x1+1);

      __m256d p1 = _mm256_mul_pd(x1,alpha);
      __m256d p2 = _mm256_mul_pd(x1,alpha_conj);

      x1 = _mm256_hsub_pd(p1,p2);

      p1 = _mm256_mul_pd(x2,alpha);
      p2 = _mm256_mul_pd(x2,alpha_conj);

      x2 = _mm256_hsub_pd(p1,p2);


      STORE_256D_ALIGNED_AS(double,_xp  ,x1);
      STORE_256D_ALIGNED_AS(double,_xp+1,x2);

      _xp += 2; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 2 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 2); i++ ){
        __m256d x = LOAD_256D_UNALIGNED_AS(double,_x+i);

        __m256d p1 = _mm256_mul_pd(x,alpha);
        __m256d p2 = _mm256_mul_pd(x,alpha_conj);

        x = _mm256_hsub_pd(p1,p2);

        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++; 
      }

      for( i = 0; i < 2 - (M % 2); i++) {
        __m256d x = _mm256_setzero_pd();
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;
      }

      _x += LDX;
    }
  }

};


// NPACK2, conj no scalar
template <>
inline void NPACKC2(const HAXX_INT M, const HAXX_INT N, 
  quaternion<double> *X, const HAXX_INT LDX, 
  quaternion<double> *Xp) {


  HAXX_INT i,j;
  quaternion<double> *_x,*_x1,*_x2;
  quaternion<double> *_xp;

  _x = X;
  _xp = Xp;

  i = M / 2;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 2;


    for( j = 0; j < N; j++ ) {

      __m256d x1 = LOAD_256D_UNALIGNED_AS(double,_x1);
      __m256d x2 = LOAD_256D_UNALIGNED_AS(double,_x1+1);

      x1 = QCONJ_256D(x1);
      x2 = QCONJ_256D(x2);


      STORE_256D_ALIGNED_AS(double,_xp  ,x1);
      STORE_256D_ALIGNED_AS(double,_xp+1,x2);

      _xp += 2; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 2 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 2); i++ ){
        __m256d x = LOAD_256D_UNALIGNED_AS(double,_x+i);
        x = QCONJ_256D(x);
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++; 
      }

      for( i = 0; i < 2 - (M % 2); i++) {
        __m256d x = _mm256_setzero_pd();
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;
      }

      _x += LDX;
    }
  }

};

// NPACK2, with conj and real scalar
template <>
inline void NPACKC2(const double ALPHA, const HAXX_INT M, 
  const HAXX_INT N, quaternion<double> *X, 
  const HAXX_INT LDX, quaternion<double> *Xp) {

  HAXX_INT i,j;
  quaternion<double> *_x,*_x1,*_x2;
  quaternion<double> *_xp;

  _x = X;
  _xp = Xp;

  const __m256d alpha = _mm256_broadcast_sd(&ALPHA);

  i = M / 2;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 2;


    for( j = 0; j < N; j++ ) {

      __m256d x1 = LOAD_256D_UNALIGNED_AS(double,_x1);
      __m256d x2 = LOAD_256D_UNALIGNED_AS(double,_x1+1);

      x1 = QCONJ_256D(x1);
      x2 = QCONJ_256D(x2);

      x1 = _mm256_mul_pd(alpha,x1);
      x2 = _mm256_mul_pd(alpha,x2);


      STORE_256D_ALIGNED_AS(double,_xp  ,x1);
      STORE_256D_ALIGNED_AS(double,_xp+1,x2);

      _xp += 2; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 2 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 2); i++ ){
        __m256d x = LOAD_256D_UNALIGNED_AS(double,_x+i);
        x = QCONJ_256D(x);
        x = _mm256_mul_pd(alpha,x);
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++; 
      }

      for( i = 0; i < 2 - (M % 2); i++) {
        __m256d x = _mm256_setzero_pd();
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;
      }

      _x += LDX;
    }
  }

};

// NPACK2, with conj and complex scalar
template <>
inline void NPACKC2(const std::complex<double> ALPHA, 
  const HAXX_INT M, const HAXX_INT N, quaternion<double> *X, 
  const HAXX_INT LDX, quaternion<double> *Xp) {

  HAXX_INT i,j;
  quaternion<double> *_x,*_x1,*_x2;
  quaternion<double> *_xp;

  _x = X;
  _xp = Xp;

  const __m256i maskConj = _mm256_set_epi64x(
                             0x8000000000000000, 0,
                             0x8000000000000000, 0 );
  
  const __m128d alphaC = LOAD_128D_UNALIGNED_AS(double,&ALPHA);
  const __m256d alpha  = SET_256D_FROM_128D(alphaC,alphaC);
  const __m256d alpha_conj = 
    _mm256_xor_pd(_mm256_permute_pd(alpha, 0x5),
                  _mm256_castsi256_pd(maskConj)
                 );

  i = M / 2;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 2;


    for( j = 0; j < N; j++ ) {

      __m256d x1 = LOAD_256D_UNALIGNED_AS(double,_x1);
      __m256d x2 = LOAD_256D_UNALIGNED_AS(double,_x1+1);

      x1 = QCONJ_256D(x1);
      x2 = QCONJ_256D(x2);

      __m256d p1 = _mm256_mul_pd(x1,alpha);
      __m256d p2 = _mm256_mul_pd(x1,alpha_conj);

      x1 = _mm256_hsub_pd(p1,p2);

      p1 = _mm256_mul_pd(x2,alpha);
      p2 = _mm256_mul_pd(x2,alpha_conj);

      x2 = _mm256_hsub_pd(p1,p2);


      STORE_256D_ALIGNED_AS(double,_xp  ,x1);
      STORE_256D_ALIGNED_AS(double,_xp+1,x2);

      _xp += 2; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 2 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 2); i++ ){
        __m256d x = LOAD_256D_UNALIGNED_AS(double,_x+i);

        x = QCONJ_256D(x);

        __m256d p1 = _mm256_mul_pd(x,alpha);
        __m256d p2 = _mm256_mul_pd(x,alpha_conj);

        x = _mm256_hsub_pd(p1,p2);

        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++; 
      }

      for( i = 0; i < 2 - (M % 2); i++) {
        __m256d x = _mm256_setzero_pd();
        STORE_256D_ALIGNED_AS(double,_xp,x);
        _xp++;
      }

      _x += LDX;
    }
  }

};




}; // namespace HAXX

