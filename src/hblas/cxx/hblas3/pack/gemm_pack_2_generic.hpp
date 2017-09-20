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
template <typename T>
inline void TPACK2(const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {

  HAXX_INT i,j;


  T *_x,*_x1,*_x2;
  T *_xp;
  

  _x = X;
  _xp = Xp;


  j = N / 2;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
 
    _x += 2*LDX;

    for( i = 0; i < M; i++ ) {
      _xp[0] = *_x1;
      _xp[1] = *_x2;

      _xp += 2;
      _x1++; _x2++;
    }
    
    j--;

  } while(j > 0);

  if( N % 2 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 2); j++) {
        _xp[0] = _x[j*LDX];
        _xp++;
      }

      for( j = 0; j < 2 - (N % 2); j++) {
        _xp[0] = 0.;
        _xp++;
      }

      _x++;
    }
  }

};

// TPACK2 with scalar
template <typename T, typename U>
inline void TPACK2(const U ALPHA, const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {

  HAXX_INT i,j;


  T *_x,*_x1,*_x2;
  T *_xp;
  

  _x = X;
  _xp = Xp;


  j = N / 2;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
 
    _x += 2*LDX;

    for( i = 0; i < M; i++ ) {
      _xp[0] = ALPHA*(*_x1);
      _xp[1] = ALPHA*(*_x2);

      _xp += 2;
      _x1++; _x2++;
    }
    
    j--;

  } while(j > 0);

  if( N % 2 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 2); j++) {
        _xp[0] = ALPHA*_x[j*LDX];
        _xp++;
      }

      for( j = 0; j < 2 - (N % 2); j++) {
        _xp[0] = 0.;
        _xp++;
      }

      _x++;
    }
  }

};


// TPACK2 with conjugation, no scalar
template <typename T>
inline void TPACKC2(const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {

  HAXX_INT i,j;


  T *_x,*_x1,*_x2;
  T *_xp;
  

  _x = X;
  _xp = Xp;


  j = N / 2;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
 
    _x += 2*LDX;

    for( i = 0; i < M; i++ ) {
      _xp[0] = SmartConj(*_x1);
      _xp[1] = SmartConj(*_x2);

      _xp += 2;
      _x1++; _x2++;
    }
    
    j--;
  } while(j > 0);

  if( N % 2 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 2); j++) {
        _xp[0] = SmartConj(_x[j*LDX]);
        _xp++;
      }

      for( j = 0; j < 2 - (N % 2); j++) {
        _xp[0] = 0.;
        _xp++;
      }

      _x++;
    }
  }

};


// TPACK2 with conjugation, with scalar
template <typename T, typename U>
inline void TPACKC2(const U ALPHA, const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {

  HAXX_INT i,j;


  T *_x,*_x1,*_x2;
  T *_xp;
  

  _x = X;
  _xp = Xp;


  j = N / 2;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
 
    _x += 2*LDX;

    for( i = 0; i < M; i++ ) {
      _xp[0] = ALPHA*SmartConj(*_x1);
      _xp[1] = ALPHA*SmartConj(*_x2);

      _xp += 2;
      _x1++; _x2++;
    }
    
    j--;
  } while(j > 0);

  if( N % 2 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 2); j++) {
        _xp[0] = ALPHA*SmartConj(_x[j*LDX]);
        _xp++;
      }

      for( j = 0; j < 2 - (N % 2); j++) {
        _xp[0] = 0.;
        _xp++;
      }

      _x++;
    }
  }

};


// NPACK2, no scalar
template <typename T>
inline void NPACK2(const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {


  HAXX_INT i,j;
  T *_x,*_x1,*_x2;
  T *_xp;

  _x = X;
  _xp = Xp;

  i = M / 2;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 2;


    for( j = 0; j < N; j++ ) {

      _xp[0] = _x1[0];
      _xp[1] = _x1[1];

      _xp += 2; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 2 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 2); i++ ){
        _xp[0] = _x[i]; _xp++; 
      }

      for( i = 0; i < 2 - (M % 2); i++) {
        _xp[0] = 0.; _xp++;
      }

      _x += LDX;
    }
  }

};

// NPACK2, with scalar
template <typename T, typename U>
inline void NPACK2(const U ALPHA, const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {


  HAXX_INT i,j;
  T *_x,*_x1,*_x2;
  T *_xp;

  _x = X;
  _xp = Xp;

  i = M / 2;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 2;


    for( j = 0; j < N; j++ ) {

      _xp[0] = ALPHA*_x1[0];
      _xp[1] = ALPHA*_x1[1];

      _xp += 2; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 2 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 2); i++ ){
        _xp[0] = ALPHA*_x[i]; _xp++; 
      }

      for( i = 0; i < 2 - (M % 2); i++) {
        _xp[0] = 0.; _xp++;
      }

      _x += LDX;
    }
  }

};


// NPACK2 with conj, no scalar
template <typename T>
inline void NPACKC2(const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {


  HAXX_INT i,j;
  T *_x,*_x1,*_x2;
  T *_xp;

  _x = X;
  _xp = Xp;

  i = M / 2;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 2;


    for( j = 0; j < N; j++ ) {

      _xp[0] = SmartConj(_x1[0]);
      _xp[1] = SmartConj(_x1[1]);

      _xp += 2; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 2 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 2); i++ ){
        _xp[0] = SmartConj(_x[i]); _xp++; 
      }

      for( i = 0; i < 2 - (M % 2); i++) {
        _xp[0] = 0.; _xp++;
      }

      _x += LDX;
    }
  }

};

// NPACK2 with conj, with scalar
template <typename T, typename U>
inline void NPACKC2(const U ALPHA, const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {


  HAXX_INT i,j;
  T *_x,*_x1,*_x2;
  T *_xp;

  _x = X;
  _xp = Xp;

  i = M / 2;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 2;


    for( j = 0; j < N; j++ ) {

      _xp[0] = ALPHA*SmartConj(_x1[0]);
      _xp[1] = ALPHA*SmartConj(_x1[1]);

      _xp += 2; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 2 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 2); i++ ){
        _xp[0] = ALPHA*SmartConj(_x[i]); _xp++; 
      }

      for( i = 0; i < 2 - (M % 2); i++) {
        _xp[0] = 0.; _xp++;
      }

      _x += LDX;
    }
  }

};


}; // namespace HAXX

