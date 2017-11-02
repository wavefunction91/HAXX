/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#include "haxx.hpp"
#include "hblas/hblas1.hpp"
#include "hblas/hblas3.hpp"

#include "util/simd.hpp"
#include "util/macro.hpp"

#include "hblas/config/types.hpp"
#include "hblas/config/hblas3/gemm.hpp"









namespace HAXX {

// Optimized GEMM
  
template <typename T, typename U, typename V>
void Kern(HAXX_INT M, HAXX_INT N, HAXX_INT K, T* __restrict__ A, U* __restrict__ B, V* __restrict__ C, 
  HAXX_INT LDC);






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

  // Packed matricies (aligned)
  _BMATF *bPack = 
    (_BMATF*)aligned_alloc(REQ_ALIGN,FixMod(KC*NC,REQ_ALIGN)*sizeof(_BMATF));
  _AMATF *aPack = 
    (_AMATF*)aligned_alloc(REQ_ALIGN,FixMod(KC*MC,REQ_ALIGN)*sizeof(_AMATF));


  // Scale C (on the left) by BETA
  HBLAS_SCALM('L','N',M,N,BETA,C,LDC,1);

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

// Turn off B Pack
#if 1

#ifdef _FACTOR_ALPHA_IN_B_PACK
      if( BTRAN )      BPACKT (nJ,nK,Bp,LDB,bPack,ALPHA);
      else if( BCT )   BPACKCT(nJ,nK,Bp,LDB,bPack,ALPHA);
      else if( BCONJ ) BPACKR (nK,nJ,Bp,LDB,bPack,ALPHA);
      else             BPACK  (nK,nJ,Bp,LDB,bPack,ALPHA);
#else
      if( BTRAN )      BPACKT (nJ,nK,Bp,LDB,bPack);
      else if( BCT )   BPACKCT(nJ,nK,Bp,LDB,bPack);
      else if( BCONJ ) BPACKR (nK,nJ,Bp,LDB,bPack);
      else             BPACK  (nK,nJ,Bp,LDB,bPack);
#endif

#endif

      for( i = 0; i < M; i += MC ) {

        nI  = std::min(M-i,MC);
        iDo = FixMod(nI,MC);

        BL1  = bPack;
        CBlk = Ci;

// Turn off A Pack
#if 1

#ifdef _FACTOR_ALPHA_IN_A_PACK
        if( ATRAN )      APACKT (nK,nI,Ai,LDA,aPack,ALPHA);
        else if( ACT )   APACKCT(nK,nI,Ai,LDA,aPack,ALPHA);
        else if( ACONJ ) APACKR (nI,nK,Ai,LDA,aPack,ALPHA);
        else             APACK  (nI,nK,Ai,LDA,aPack,ALPHA);
#else
        if( ATRAN )      APACKT (nK,nI,Ai,LDA,aPack);
        else if( ACT )   APACKCT(nK,nI,Ai,LDA,aPack);
        else if( ACONJ ) APACKR (nI,nK,Ai,LDA,aPack);
        else             APACK  (nI,nK,Ai,LDA,aPack);
#endif

#endif


// Turn off Kernel evaluation
#if 1
        for( jj = 0; jj < nJ; jj += NR ) {
          nJJJ = std::min(NR,nJ-jj);

          smallC = CBlk;
          smallA = aPack;

          for( ii = 0; ii < nI; ii += MR ) {
            nIII = std::min(MR,nI-ii);

            // Perform kernel operation
            Kern(nIII,nJJJ,nK,smallA,BL1,smallC,LDC);
     
            smallC += MR;
            smallA += MR*nK;
          }

          BL1  += NR*nK;
          CBlk += NR*LDC;
        }
#endif

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

