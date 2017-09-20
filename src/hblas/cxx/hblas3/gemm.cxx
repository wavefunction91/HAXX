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

#include <algorithm>

#define FixMod(X,N) (( (X) % (N) ) ? (X) + (N) - ((X) % (N)) : (X))

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
#define NR 2

#if MR != 4 && MR != 2
  #error MR must be 4 or 2
#endif

#if NR != 4 && NR != 2
  #error NR must be 4 or 2
#endif

#include "gemm_pack_4.hpp"
#include "gemm_pack_2.hpp"

#define _FACTOR_ALPHA_IN_A_PACK
//#define _FACTOR_ALPHA_IN_B_PACK


namespace HAXX {

// Optimized GEMM
  

template <typename T, typename U, typename V, typename _BetaF>
inline void Kern(_BetaF BETA, HAXX_INT M, HAXX_INT N, HAXX_INT K, T* __restrict__ A, U* __restrict__ B, V* __restrict__ C, 
  HAXX_INT LDC);

template <typename T, typename U, typename V, typename _BetaF>
inline void Kern(_BetaF BETA, HAXX_INT M, HAXX_INT N, HAXX_INT K, T*  A, U*  B, V*  C, 
  HAXX_INT LDC) {

  for(auto j = 0; j < N; j++)
  for(auto i = 0; i < M; i++)
    C[i + j*LDC] *= BETA;

  for(auto k = 0; k < K; k++) {

    for(auto j = 0; j < N; j++)
    for(auto i = 0; i < M; i++)
      C[i + j*LDC] += A[i] * B[j];

    A += MR;
    B += NR;

  }

}


#if 0
#if MR == 4 && NR == 2
template <>
inline void Kern(double BETA, HAXX_INT M, HAXX_INT N, HAXX_INT K, 
  quaternion<double>* restrict A, quaternion<double>* restrict B, quaternion<double>* restrict C, 
  HAXX_INT LDC) {

  // Load 4x2 block of C
  __m256d c00 = LOADD_UNALIGNED_AS(double,C  );
  __m256d c10 = LOADD_UNALIGNED_AS(double,C+1);
  __m256d c20 = LOADD_UNALIGNED_AS(double,C+2);
  __m256d c30 = LOADD_UNALIGNED_AS(double,C+3);

  __m256d c01 = LOADD_UNALIGNED_AS(double,C+LDC  );
  __m256d c11 = LOADD_UNALIGNED_AS(double,C+LDC+1);
  __m256d c21 = LOADD_UNALIGNED_AS(double,C+LDC+2);
  __m256d c31 = LOADD_UNALIGNED_AS(double,C+LDC+3);

}
#endif
#endif






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

  std::cout << "HEREX\n";

  // Packed matricies (aligned)
  _BMATF *bPack = 
    (_BMATF*)aligned_alloc(REQ_ALIGN,FixMod(KC*NC,REQ_ALIGN)*sizeof(_BMATF));
  _AMATF *aPack = 
    (_AMATF*)aligned_alloc(REQ_ALIGN,FixMod(KC*MC,REQ_ALIGN)*sizeof(_AMATF));


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
    jDo = FixMod(nJ,NR);

    Bp = Bj;
    Ap = A;

    for( k = 0; k < K; k += KC ) {

      nK = std::min(K-k,KC);

      Ai = Ap;
      Ci = Cj;

#ifdef _FACTOR_ALPHA_IN_B_PACK
  #if NR == 4
      if( BTRAN )      NPACK4 (ALPHA,nJ,nK,Bp,LDB,bPack);
      else if( BCT )   NPACKC4(ALPHA,nJ,nK,Bp,LDB,bPack);
      else if( BCONJ ) TPACKC4(ALPHA,nK,nJ,Bp,LDB,bPack);
      else             TPACK4 (ALPHA,nK,nJ,Bp,LDB,bPack);
  #elif NR == 2
      if( BTRAN )      NPACK2 (ALPHA,nJ,nK,Bp,LDB,bPack);
      else if( BCT )   NPACKC2(ALPHA,nJ,nK,Bp,LDB,bPack);
      else if( BCONJ ) TPACKC2(ALPHA,nK,nJ,Bp,LDB,bPack);
      else             TPACK2 (ALPHA,nK,nJ,Bp,LDB,bPack);
  #endif
#else
  #if NR == 4
      if( BTRAN )      NPACK4 (nJ,nK,Bp,LDB,bPack);
      else if( BCT )   NPACKC4(nJ,nK,Bp,LDB,bPack);
      else if( BCONJ ) TPACKC4(nK,nJ,Bp,LDB,bPack);
      else             TPACK4 (nK,nJ,Bp,LDB,bPack);
  #elif NR == 2
      if( BTRAN )      NPACK2 (nJ,nK,Bp,LDB,bPack);
      else if( BCT )   NPACKC2(nJ,nK,Bp,LDB,bPack);
      else if( BCONJ ) TPACKC2(nK,nJ,Bp,LDB,bPack);
      else             TPACK2 (nK,nJ,Bp,LDB,bPack);
  #endif
#endif

      for( i = 0; i < M; i += MC ) {

        nI  = std::min(M-i,MC);
        iDo = FixMod(nI,MC);

        BL1  = bPack;
        CBlk = Ci;

#ifdef _FACTOR_ALPHA_IN_A_PACK
  #if MR == 4
        if( ATRAN )      TPACK4 (ALPHA,nK,nI,Ai,LDA,aPack);
        else if( ACT )   TPACKC4(ALPHA,nK,nI,Ai,LDA,aPack);
        else if( ACONJ ) NPACKC4(ALPHA,nI,nK,Ai,LDA,aPack);
        else             NPACK4 (ALPHA,nI,nK,Ai,LDA,aPack);
  #elif MR == 2
        if( ATRAN )      TPACK2 (ALPHA,nK,nI,Ai,LDA,aPack);
        else if( ACT )   TPACKC2(ALPHA,nK,nI,Ai,LDA,aPack);
        else if( ACONJ ) NPACKC2(ALPHA,nI,nK,Ai,LDA,aPack);
        else             NPACK2 (ALPHA,nI,nK,Ai,LDA,aPack);
  #endif
#else
  #if MR == 4
        if( ATRAN )      TPACK4 (nK,nI,Ai,LDA,aPack);
        else if( ACT )   TPACKC4(nK,nI,Ai,LDA,aPack);
        else if( ACONJ ) NPACKC4(nI,nK,Ai,LDA,aPack);
        else             NPACK4 (nI,nK,Ai,LDA,aPack);
  #elif MR == 2
        if( ATRAN )      TPACK2 (nK,nI,Ai,LDA,aPack);
        else if( ACT )   TPACKC2(nK,nI,Ai,LDA,aPack);
        else if( ACONJ ) NPACKC2(nI,nK,Ai,LDA,aPack);
        else             NPACK2 (nI,nK,Ai,LDA,aPack);
  #endif
  
  #ifndef _FACTOR_ALPHA_IN_B_PACK
        std::transform(aPack,&aPack[nK*iDo],aPack,[&](_AMATF x){ return ALPHA*x;});
  #endif
#endif


        for( jj = 0; jj < jDo; jj += NR ) {

          nJJJ = std::min(NR,nJ-jj);

          smallC = CBlk;
          smallA = aPack;

          for( ii = 0; ii < iDo; ii += MR ) {
            nIII = std::min(MR,nI-ii);

            // Perform kernel operation
            Kern(BETA,nIII,nJJJ,nK,smallA,BL1,smallC,LDC);
     
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

  // Free packed matricies
  free(bPack); free(aPack);

}; // HBLAS_GEMM


}; // namespace HAXX

