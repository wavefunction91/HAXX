/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#include "haxx.hpp"
#include "hblas/hblas3_def.hpp"

#include "util/simd.hpp"

#include <vector>
#include <algorithm>

// Determine type of scaling parameter ALPHA
#if ALPHAF == DOUBLE
  #define _ALPHAF double
#elif ALPHAF == DCOMPLEX
  #define _ALPHAF std::complex<double>
#elif ALPHAF == DQUATERNION
  #define _ALPHAF quaternion<double>
#else
  #error GEMM Only Supports 64-bit floats
#endif

// Determine type of scaling parameter BETA
#if BETAF == DOUBLE
  #define _BETAF double
#elif BETAF == DCOMPLEX
  #define _BETAF std::complex<double>
#elif BETAF == DQUATERNION
  #define _BETAF quaternion<double>
#else
  #error GEMM Only Supports 64-bit floats
#endif

#define _AMATF quaternion<double>
#define _BMATF quaternion<double>


#define MC 128
#define NC 512
#define KC 128
#define MR 4
#define NR 4


namespace HAXX {

// Optimized GEMM
  

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

template <typename T, typename U, typename V>
void Kern(HAXX_INT K, T* A, U* B, V* C, HAXX_INT LDC) {

  for(auto k = 0; k < K; k++) {

    for(auto j = 0; j < NR; j++)
    for(auto i = 0; i < MR; i++)
      C[i + j*LDC] += A[i] * B[j];

    A += MR;
    B += NR;

  }

}


template<>
void HBLAS_GEMM(const char TRANSA, const char TRANSB, const HAXX_INT M, 
  const HAXX_INT N, const HAXX_INT K, const _ALPHAF ALPHA, _AMATF * const A, 
  const HAXX_INT LDA, _BMATF * const B, const HAXX_INT LDB, 
  const _BETAF BETA, quaternion<double> * const C, const HAXX_INT LDC) {


  const bool ATRAN = TRANSA == 'T';
  const bool ACT   = TRANSA == 'C';
  const bool ACONJ = TRANSA == 'R';

  const bool BTRAN = TRANSB == 'T';
  const bool BCT   = TRANSB == 'C';
  const bool BCONJ = TRANSB == 'R';


  // Packed matricies
  std::vector<_BMATF> bPack(KC*NC);
  std::vector<_AMATF> aPack(MC*KC);
  std::vector<quaternion<double>> cPack(NR*MR);

  // Counter vars
  HAXX_INT j, nJ, jDo;
  HAXX_INT i, nI, iDo;
  HAXX_INT k, nK;
  HAXX_INT ii, jj, iii, jjj, nIII, nJJJ;

  // Initial N-cut panels of C and B
  quaternion<double> *Cj = C;
  _BMATF             *Bj = B;


  // Panel / block pointers
  _BMATF *Bp, *BL1;
  _AMATF *Ap, *Ai, *smallA;
  quaternion<double> *Ci, *CBlk, *smallC;

  for( j = 0; j < N; j += NC ) {

    nJ  = std::min(N-j,NC);
    jDo = ( nJ % NR ) ? nJ + NR - (nJ % NR) : nJ;


    Bp = Bj;
    Ap = A;

    for( k = 0; k < K; k += KC ) {

      nK = std::min(K-k,KC);

      Ai = Ap;
      Ci = Cj;

      if( BTRAN )      NPACK4 (nJ,nK,Bp,LDB,&bPack[0]);
      else if( BCT )   NPACKC4(nJ,nK,Bp,LDB,&bPack[0]);
      else if( BCONJ ) TPACKC4(nK,nJ,Bp,LDB,&bPack[0]);
      else             TPACK4 (nK,nJ,Bp,LDB,&bPack[0]);

      for( i = 0; i < M; i += MC ) {

        nI = std::min(M-i,MC);
        iDo = ( nI % MR ) ? nI + MR - (nI % MR): nI;

        BL1  = &bPack[0];
        CBlk = Ci;

        if( ATRAN )      TPACK4 (nK,nI,Ai,LDA,&aPack[0]);
        else if( ACT )   TPACKC4(nK,nI,Ai,LDA,&aPack[0]);
        else if( ACONJ ) NPACKC4(nI,nK,Ai,LDA,&aPack[0]);
        else             NPACK4 (nI,nK,Ai,LDA,&aPack[0]);

        // FIXME: terrible idea...
        std::transform(&aPack[0],&aPack[nK*iDo],&aPack[0],[&](_AMATF x){ return ALPHA*x;});

        for( jj = 0; jj < jDo; jj += NR ) {

          nJJJ = std::min(NR,nJ-jj);

          smallC = CBlk;
          smallA = &aPack[0];

          for( ii = 0; ii < iDo; ii += MR ) {
            nIII = std::min(MR,nI-ii);

            for( jjj = 0; jjj < nJJJ; jjj++ )
            for( iii = 0; iii < nIII; iii++ ) 
              cPack[iii + jjj*MR] = BETA*smallC[iii +jjj*LDC];

            // Perform kernel operation
            Kern(nK,smallA,BL1,&cPack[0],MR);

            for( jjj = 0; jjj < nJJJ; jjj++ )
            for( iii = 0; iii < nIII; iii++ ) 
              smallC[iii +jjj*LDC] = cPack[iii + jjj*MR];

            smallC += MR;
            smallA += MR*nK;

          }

          BL1  += NR*nK;

          CBlk += NR*LDC;
        }

        if( ATRAN or ACT ) Ai += nI*LDA;
        else               Ai += nI; 

        Ci += nI;
      }

      if( ATRAN or ACT ) Ap += nK;
      else               Ap += nK*LDA;

      if( BTRAN or BCT ) Bp += nK*LDB;
      else               Bp += nK;
    }

    Cj += nJ*LDC;

    if( BTRAN or BCT ) Bj += nJ;
    else               Bj += nJ*LDB;
  }


}; // HBLAS_GEMM


}; // namespace HAXX

