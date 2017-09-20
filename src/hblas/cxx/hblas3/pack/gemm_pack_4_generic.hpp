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


// TPACK4, no scalar
template <typename T>
inline void TPACK4(const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {

  HAXX_INT i,j;


  T *_x,*_x1,*_x2,*_x3,*_x4;
  T *_xp;
  

  _x = X;
  _xp = Xp;


  j = N / 4;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
    _x3 = _x2 + LDX;
    _x4 = _x3 + LDX;
 
    _x += 4*LDX;

    for( i = 0; i < M; i++ ) {
      _xp[0] = *_x1;
      _xp[1] = *_x2;
      _xp[2] = *_x3;
      _xp[3] = *_x4;

      _xp += 4;
      _x1++; _x2++; _x3++; _x4++;
    }
    
    j--;

  } while(j > 0);

  if( N % 4 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 4); j++) {
        _xp[0] = _x[j*LDX];
        _xp++;
      }

      for( j = 0; j < 4 - (N % 4); j++) {
        _xp[0] = 0.;
        _xp++;
      }

      _x++;
    }
  }

};

// TPACK4 with scalar
template <typename T, typename U>
inline void TPACK4(const U ALPHA, const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {

  HAXX_INT i,j;


  T *_x,*_x1,*_x2,*_x3,*_x4;
  T *_xp;
  

  _x = X;
  _xp = Xp;


  j = N / 4;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
    _x3 = _x2 + LDX;
    _x4 = _x3 + LDX;
 
    _x += 4*LDX;

    for( i = 0; i < M; i++ ) {
      _xp[0] = ALPHA*(*_x1);
      _xp[1] = ALPHA*(*_x2);
      _xp[2] = ALPHA*(*_x3);
      _xp[3] = ALPHA*(*_x4);

      _xp += 4;
      _x1++; _x2++; _x3++; _x4++;
    }
    
    j--;

  } while(j > 0);

  if( N % 4 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 4); j++) {
        _xp[0] = ALPHA*_x[j*LDX];
        _xp++;
      }

      for( j = 0; j < 4 - (N % 4); j++) {
        _xp[0] = 0.;
        _xp++;
      }

      _x++;
    }
  }

};


// TPACK4 with conj, no scalar
template <typename T>
inline void TPACKC4(const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {

  HAXX_INT i,j;


  T *_x,*_x1,*_x2,*_x3,*_x4;
  T *_xp;
  

  _x = X;
  _xp = Xp;


  j = N / 4;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
    _x3 = _x2 + LDX;
    _x4 = _x3 + LDX;
 
    _x += 4*LDX;

    for( i = 0; i < M; i++ ) {
      _xp[0] = SmartConj(*_x1);
      _xp[1] = SmartConj(*_x2);
      _xp[2] = SmartConj(*_x3);
      _xp[3] = SmartConj(*_x4);

      _xp += 4;
      _x1++; _x2++; _x3++; _x4++;
    }
    
    j--;
  } while(j > 0);

  if( N % 4 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 4); j++) {
        _xp[0] = SmartConj(_x[j*LDX]);
        _xp++;
      }

      for( j = 0; j < 4 - (N % 4); j++) {
        _xp[0] = 0.;
        _xp++;
      }

      _x++;
    }
  }

};


// TPACK4 with conj, with scalar
template <typename T, typename U>
inline void TPACKC4(const U ALPHA, const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {

  HAXX_INT i,j;


  T *_x,*_x1,*_x2,*_x3,*_x4;
  T *_xp;
  

  _x = X;
  _xp = Xp;


  j = N / 4;
  if( j > 0 )
  do {
    _x1 = _x;
    _x2 = _x1 + LDX;
    _x3 = _x2 + LDX;
    _x4 = _x3 + LDX;
 
    _x += 4*LDX;

    for( i = 0; i < M; i++ ) {
      _xp[0] = ALPHA*SmartConj(*_x1);
      _xp[1] = ALPHA*SmartConj(*_x2);
      _xp[2] = ALPHA*SmartConj(*_x3);
      _xp[3] = ALPHA*SmartConj(*_x4);

      _xp += 4;
      _x1++; _x2++; _x3++; _x4++;
    }
    
    j--;
  } while(j > 0);

  if( N % 4 ) {
    for( i = 0; i < M; i++ ) {

      for( j = 0; j < (N % 4); j++) {
        _xp[0] = ALPHA*SmartConj(_x[j*LDX]);
        _xp++;
      }

      for( j = 0; j < 4 - (N % 4); j++) {
        _xp[0] = 0.;
        _xp++;
      }

      _x++;
    }
  }

};


// NPACK4, no scalar
template <typename T>
inline void NPACK4(const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {


  HAXX_INT i,j;
  T *_x,*_x1,*_x2,*_x3,*_x4;
  T *_xp;

  _x = X;
  _xp = Xp;

  i = M / 4;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 4;


    for( j = 0; j < N; j++ ) {

      _xp[0] = _x1[0];
      _xp[1] = _x1[1];
      _xp[2] = _x1[2];
      _xp[3] = _x1[3];

      _xp += 4; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 4 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 4); i++ ){
        _xp[0] = _x[i]; _xp++; 
      }

      for( i = 0; i < 4 - (M % 4); i++) {
        _xp[0] = 0.; _xp++;
      }

      _x += LDX;
    }
  }

};


// NPACK4, with scalar
template <typename T,typename U>
inline void NPACK4(const U ALPHA, const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {


  HAXX_INT i,j;
  T *_x,*_x1,*_x2,*_x3,*_x4;
  T *_xp;

  _x = X;
  _xp = Xp;

  i = M / 4;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 4;


    for( j = 0; j < N; j++ ) {

      _xp[0] = ALPHA*_x1[0];
      _xp[1] = ALPHA*_x1[1];
      _xp[2] = ALPHA*_x1[2];
      _xp[3] = ALPHA*_x1[3];

      _xp += 4; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 4 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 4); i++ ){
        _xp[0] = ALPHA*_x[i]; _xp++; 
      }

      for( i = 0; i < 4 - (M % 4); i++) {
        _xp[0] = 0.; _xp++;
      }

      _x += LDX;
    }
  }

};


// NPACK 4 with conj, no scalar
template <typename T>
inline void NPACKC4(const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {


  HAXX_INT i,j;
  T *_x,*_x1,*_x2,*_x3,*_x4;
  T *_xp;

  _x = X;
  _xp = Xp;

  i = M / 4;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 4;


    for( j = 0; j < N; j++ ) {

      _xp[0] = SmartConj(_x1[0]);
      _xp[1] = SmartConj(_x1[1]);
      _xp[2] = SmartConj(_x1[2]);
      _xp[3] = SmartConj(_x1[3]);

      _xp += 4; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 4 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 4); i++ ){
        _xp[0] = SmartConj(_x[i]); _xp++; 
      }

      for( i = 0; i < 4 - (M % 4); i++) {
        _xp[0] = 0.; _xp++;
      }

      _x += LDX;
    }
  }

};

// NPACK 4 with conj, with scalar
template <typename T, typename U>
inline void NPACKC4(const U ALPHA, const HAXX_INT M, const HAXX_INT N, T *X, 
  const HAXX_INT LDX, T *Xp) {


  HAXX_INT i,j;
  T *_x,*_x1,*_x2,*_x3,*_x4;
  T *_xp;

  _x = X;
  _xp = Xp;

  i = M / 4;
  if( i > 0 )
  do {

    _x1 = _x;
    _x += 4;


    for( j = 0; j < N; j++ ) {

      _xp[0] = ALPHA*SmartConj(_x1[0]);
      _xp[1] = ALPHA*SmartConj(_x1[1]);
      _xp[2] = ALPHA*SmartConj(_x1[2]);
      _xp[3] = ALPHA*SmartConj(_x1[3]);

      _xp += 4; _x1 += LDX;
    }

    i--;
  } while( i > 0 );

  if( M % 4 ) {
    for( j = 0; j < N; j++ ) {
      for( i = 0; i < (M % 4); i++ ){
        _xp[0] = ALPHA*SmartConj(_x[i]); _xp++; 
      }

      for( i = 0; i < 4 - (M % 4); i++) {
        _xp[0] = 0.; _xp++;
      }

      _x += LDX;
    }
  }

};


}; // namespace HAXX
