/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#include "haxx.hpp"
#include "hblas/hblas1/hblas_scalm_impl.hpp"
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


#define MC 64
#define NC 1024
#define KC 64
#define MR 2
#define NR 2

#if MR != 4 && MR != 2
  #error MR must be 4 or 2
#endif

#if NR != 4 && NR != 2
  #error NR must be 4 or 2
#endif


#define _FACTOR_TRANSPOSE_INTO_PACK

#include "gemm_pack_4.hpp"
#include "gemm_pack_2.hpp"

#define _FACTOR_ALPHA_IN_A_PACK
//#define _FACTOR_ALPHA_IN_B_PACK


namespace HAXX {

// Optimized GEMM
  


template <typename T, typename U, typename V>
inline void Kern(HAXX_INT M, HAXX_INT N, HAXX_INT K, T* __restrict__ A, U* __restrict__ B, V* __restrict__ C, 
  HAXX_INT LDC) {


  for(auto k = 0; k < K; k++) {

    for(auto j = 0; j < N; j++)
    for(auto i = 0; i < M; i++)
      C[i + j*LDC] += A[i] * B[j];

    A += MR;
    B += NR;

  }

}


#if MR == 2 && NR == 2
//#if 0
inline void Kern(HAXX_INT M, HAXX_INT N, HAXX_INT K, 
  quaternion<double>* __restrict__ A, quaternion<double>* __restrict__ B, 
  quaternion<double>* __restrict__ C, HAXX_INT LDC) {

  __m256d t1,t2,t3,t4;

  // Load C
  __m256d c00 = LOAD_256D_UNALIGNED_AS(double,C      );
  __m256d c10 = LOAD_256D_UNALIGNED_AS(double,C+1    );
  __m256d c01 = LOAD_256D_UNALIGNED_AS(double,C+LDC  );
  __m256d c11 = LOAD_256D_UNALIGNED_AS(double,C+LDC+1);

  _MM_TRANSPOSE_4x4_PD(c00,c01,c10,c11,t1,t2,t3,t4);


  quaternion<double> *locB = B, *locA = A;

  HAXX_INT k = K;

  if( k > 0 )
  do {

    // Load A 
    __m256d a00 = LOAD_256D_ALIGNED_AS(double,locA);
    __m256d a10 = LOAD_256D_ALIGNED_AS(double,locA+1);
    locA += 2;
    

    // Load B 
    __m256d b00 = LOAD_256D_ALIGNED_AS(double,locB);
    __m256d b10 = LOAD_256D_ALIGNED_AS(double,locB+1);
    locB += 2;


#if 1

    __m256d a00c = a00;
    __m256d a10c = a10;
    _MM_TRANSPOSE_4x4_PD(a00,a00c,a10,a10c,t1,t2,t3,t4);

    __m256d b00c = b00;
    __m256d b10c = b10;
    _MM_TRANSPOSE_4x4_PD(b00,b10,b00c,b10c,t1,t2,t3,t4);

#else
  
    //__m256d a_SISI = _mm256_permute2f128_pd(a00,a10, 0x20);
    //__m256d a_JKJK = _mm256_permute2f128_pd(a00,a10, 0x31);
  
    __m256d a_SISI = _mm256_undefined_pd();
    __m256d a_JKJK = _mm256_undefined_pd();

    a00    = _mm256_permute_pd(a_SISI, 0x0);
    a_SISI = _mm256_permute_pd(a_SISI, 0xF);
  
    a10    = _mm256_permute_pd(a_JKJK, 0x0);
    a_JKJK = _mm256_permute_pd(a_JKJK, 0xF);
  
    __m256d &a00c = a_SISI;
    __m256d &a10c = a_JKJK;

    __m256d b00c = b00;
    __m256d b10c = b10;
    _MM_TRANSPOSE_4x4_PD(b00,b10,b00c,b10c,t1,t2,t3,t4);

#endif

    INC_MULD4Q_NN(a00,a00c,a10,a10c,b00,b10,b00c,b10c,c00,c01,c10,c11);
  

    k--;
  } while( k > 0 );

  _MM_TRANSPOSE_4x4_PD(c00,c01,c10,c11,t1,t2,t3,t4);

  STORE_256D_UNALIGNED_AS(double,C      ,c00);
  STORE_256D_UNALIGNED_AS(double,C+1    ,c10);
  STORE_256D_UNALIGNED_AS(double,C+LDC  ,c01);
  STORE_256D_UNALIGNED_AS(double,C+LDC+1,c11);
  

}
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

  // Packed matricies (aligned)
  _BMATF *bPack = 
    (_BMATF*)aligned_alloc(REQ_ALIGN,FixMod(KC*NC,REQ_ALIGN)*sizeof(_BMATF));
  _AMATF *aPack = 
    (_AMATF*)aligned_alloc(REQ_ALIGN,FixMod(KC*MC,REQ_ALIGN)*sizeof(_AMATF));

//std::cout << "BPack " << FixMod(KC*NC,REQ_ALIGN)*sizeof(_BMATF) << std::endl;
//std::cout << "APack " << FixMod(KC*MC,REQ_ALIGN)*sizeof(_AMATF) << std::endl;


  // Scale C by BETA
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

#endif

      for( i = 0; i < M; i += MC ) {

        nI  = std::min(M-i,MC);
        iDo = FixMod(nI,MC);

        BL1  = bPack;
        CBlk = Ci;

// Turn off A Pack
#if 1

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

#endif


// Turn off Kernel evaluation
#if 1
        for( jj = 0; jj < nJ; jj += NR ) {

          //nJJJ = std::min(NR,jj-nJ);
          nJJJ = NR;

          smallC = CBlk;
          smallA = aPack;

          for( ii = 0; ii < nI; ii += MR ) {
            //nIII = std::min(MR,ii-nI);
            nIII = MR;

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

