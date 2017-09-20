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
#define NR 2

#include "gemm_pack_4.hpp"
#include "gemm_pack_2.hpp"

namespace HAXX {

// Optimized GEMM
  


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

  std::cout << "HEREX\n";

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

#if NR == 4
      if( BTRAN )      NPACK4 (nJ,nK,Bp,LDB,&bPack[0]);
      else if( BCT )   NPACKC4(nJ,nK,Bp,LDB,&bPack[0]);
      else if( BCONJ ) TPACKC4(nK,nJ,Bp,LDB,&bPack[0]);
      else             TPACK4 (nK,nJ,Bp,LDB,&bPack[0]);
#elif NR == 2
      if( BTRAN )      NPACK2 (nJ,nK,Bp,LDB,&bPack[0]);
      else if( BCT )   NPACKC2(nJ,nK,Bp,LDB,&bPack[0]);
      else if( BCONJ ) TPACKC2(nK,nJ,Bp,LDB,&bPack[0]);
      else             TPACK2 (nK,nJ,Bp,LDB,&bPack[0]);
#else
  #error Something has gone horribly wrong!
#endif

      for( i = 0; i < M; i += MC ) {

        nI = std::min(M-i,MC);
        iDo = ( nI % MR ) ? nI + MR - (nI % MR): nI;

        BL1  = &bPack[0];
        CBlk = Ci;

#if MR == 4
        if( ATRAN )      TPACK4 (nK,nI,Ai,LDA,&aPack[0]);
        else if( ACT )   TPACKC4(nK,nI,Ai,LDA,&aPack[0]);
        else if( ACONJ ) NPACKC4(nI,nK,Ai,LDA,&aPack[0]);
        else             NPACK4 (nI,nK,Ai,LDA,&aPack[0]);
#elif MR == 2
        if( ATRAN )      TPACK2 (nK,nI,Ai,LDA,&aPack[0]);
        else if( ACT )   TPACKC2(nK,nI,Ai,LDA,&aPack[0]);
        else if( ACONJ ) NPACKC2(nI,nK,Ai,LDA,&aPack[0]);
        else             NPACK2 (nI,nK,Ai,LDA,&aPack[0]);
#else
  #error Something has gone horribly wrong!
#endif

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

