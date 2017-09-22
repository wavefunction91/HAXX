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
#define MR 2
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


template <typename _BetaF>
inline void scaleC(_BetaF &BETA, __m256d &c1, __m256d &c2, __m256d &c3, __m256d &c4);

template <>
inline void scaleC(double &BETA, __m256d &c1, __m256d &c2, __m256d &c3, __m256d &c4) {

  const __m256d  beta = _mm256_broadcast_sd(&BETA);
  c1 = _mm256_mul_pd(beta,c1);
  c2 = _mm256_mul_pd(beta,c2);
  c3 = _mm256_mul_pd(beta,c3);
  c4 = _mm256_mul_pd(beta,c4);

};

template <>
inline void scaleC(std::complex<double> &BETA, __m256d &c1, __m256d &c2, __m256d &c3, __m256d &c4) {

  const __m256i maskConj = _mm256_set_epi64x(
                             0x8000000000000000, 0,
                             0x8000000000000000, 0 );

  __m128d beta = _mm_loadu_pd(reinterpret_cast<double*>(&BETA));

  // XXX: this is not correct for right multiply
  __m256d betaV  = SET_256D_FROM_128D(beta,beta);
  __m256d betaVC = _mm256_xor_pd( 
                     _mm256_permute_pd(betaV,0x5),
                     _mm256_castsi256_pd(maskConj)
                   );

  c1 = _mm256_hsub_pd(_mm256_mul_pd(c1,betaV),_mm256_mul_pd(c1,betaVC));
  c2 = _mm256_hsub_pd(_mm256_mul_pd(c2,betaV),_mm256_mul_pd(c2,betaVC));
  c3 = _mm256_hsub_pd(_mm256_mul_pd(c3,betaV),_mm256_mul_pd(c3,betaVC));
  c4 = _mm256_hsub_pd(_mm256_mul_pd(c4,betaV),_mm256_mul_pd(c4,betaVC));


};



#if MR == 2 && NR == 2
template <typename _BetaF>
inline void Kern(_BetaF BETA, HAXX_INT M, HAXX_INT N, HAXX_INT K, 
  quaternion<double>* __restrict__ A, quaternion<double>* __restrict__ B, 
  quaternion<double>* __restrict__ C, HAXX_INT LDC) {

  // Scratch registers
  // ymm12->15
  __m256d t1,t2,t3,t4;
  // xmm12->15
  __m128d x1,x2,x3,x4;


  // Load C
  //   c00 (ymm0) = [ C00(S) C00(I) C00(J) C00(K) ]
  //   c01 (ymm1) = [ C01(S) C01(I) C01(J) C01(K) ]
  //   c10 (ymm2) = [ C10(S) C10(I) C10(J) C10(K) ]
  //   c11 (ymm3) = [ C11(S) C11(I) C11(J) C11(K) ]
  __m256d c00 = LOAD_256D_UNALIGNED_AS(double,C      );
  __m256d c10 = LOAD_256D_UNALIGNED_AS(double,C+1    );
  __m256d c01 = LOAD_256D_UNALIGNED_AS(double,C+LDC  );
  __m256d c11 = LOAD_256D_UNALIGNED_AS(double,C+LDC+1);

  // Scale C by Beta
  scaleC(BETA,c00,c10,c01,c11);


  quaternion<double> *locB = B, *locA = A;


  HAXX_INT k = K / 2;

  if( k > 0 )
  do {

    // Load A (4 new registers, ymm4->7 for A)
    __m256d a00 = LOAD_256D_ALIGNED_AS(double,locA);
    __m256d a10 = LOAD_256D_ALIGNED_AS(double,locA+1);
    __m256d a01 = LOAD_256D_ALIGNED_AS(double,locA+2);
    __m256d a11 = LOAD_256D_ALIGNED_AS(double,locA+3);
    locA += 4;
    

    // Load B (4 new registers -> ymm8->11 for B)
    __m256d b00 = LOAD_256D_ALIGNED_AS(double,locB);
    __m256d b10 = LOAD_256D_ALIGNED_AS(double,locB+1);
    __m256d b01 = LOAD_256D_ALIGNED_AS(double,locB+2);
    __m256d b11 = LOAD_256D_ALIGNED_AS(double,locB+3);
    locB += 4;

    // Transpose A (8 registers, ymm4->7 for A, ymm12->15 for T)
    //   a00 (ymm4) = [ A00(S) A01(S) A10(S) A11(S) ]
    //   a01 (ymm5) = [ A00(I) A01(I) A10(I) A11(I) ]
    //   a10 (ymm6) = [ A00(J) A01(J) A10(J) A11(J) ]
    //   a11 (ymm7) = [ A00(K) A01(K) A10(K) A11(K) ]
    _MM_TRANSPOSE_4x4_PD(a00,a01,a10,a11,t1,t2,t3,t4);

    // Transpose B (8 registers, ymm8->11 for A, ymm12->15 for T)
    //   b00 (ymm8)  = [ B00(S) B01(S) B10(S) B11(S) ]
    //   b01 (ymm9)  = [ B00(I) B01(I) B10(I) B11(I) ]
    //   b10 (ymm10) = [ B00(J) B01(J) B10(J) B11(J) ]
    //   b11 (ymm11) = [ B00(K) B01(K) B10(K) B11(K) ]
    _MM_TRANSPOSE_4x4_PD(b00,b01,b10,b11,t1,t2,t3,t4);
    


    // Let PIJKL = AIJ * BKL
    //  s.t.
    //    C00 = P0000 + P0101
    //    C10 = P1000 + P1101
    //    C01 = P0010 + P0111
    //    C11 = P1010 + P1111



    // COMPUTE THE DIAGONAL CONTRIBUTIONS TO THE PRODUCT



    // Compute product in T (ymm12-ymm15)
    // t1 (ymm12) = [ P0000(S) P0101(S) P1010(S) P1111(S) ]
    // t2 (ymm13) = [ P0000(I) P0101(I) P1010(I) P1111(I) ]
    // t3 (ymm14) = [ P0000(J) P0101(J) P1010(J) P1111(J) ]
    // t4 (ymm15) = [ P0000(K) P0101(K) P1010(K) P1111(K) ]
    // Multiplication requires no register stratch?
    MULD4Q_NN(a00,a01,a10,a11,b00,b01,b10,b11,t1,t2,t3,t4);


    // t1 (ymm12) =                t3 (ymm14) = 
    // [                           [
    //   P0000(S) + P0101(S)         P0000(J) + P0101(J)
    //   P0000(I) + P0101(I)         P0000(K) + P0101(K)
    //   P1010(S) + P1111(S)         P1010(J) + P1111(J)
    //   P1010(I) + P1111(I)         P1010(K) + P1111(K)
    // ]                           ]
    t1 = _mm256_hadd_pd(t1,t2);
    t3 = _mm256_hadd_pd(t3,t4);


    // XXX: Is there a better way to set the high bits of one
    // register to the low bits of another?


    // x2 (xmm13) =              x4 (xmm15) =
    // [                         [
    //   P0000(S) + P0101(S)       P0000(J) + P0101(J)
    //   P0000(I) + P0101(I)       P0000(K) + P0101(K)
    // ]                         ]
    x2 = GET_LO_128D_256D(t1);
    x4 = GET_LO_128D_256D(t3);


    // t2 (ymm13) =
    // [
    //   P0000(S) + P0101(S)
    //   P0000(I) + P0101(I)
    //   P0000(J) + P0101(J)
    //   P0000(K) + P0101(K)
    // ]
    t2 = SET_256D_FROM_128D(x2,x4); 
    
    // C00 (ymm0) += t2 (ymm13)
    c00 = _mm256_add_pd(c00,t2);


    // x2 (xmm13) =             x4 (xmm15) =
    // [                        [
    //   P1010(S) + P1111(S)      P1010(J) + P1111(J)
    //   P1010(I) + P1111(I)      P1010(K) + P1111(K)
    // ]                        ]
    x2 = GET_HI_128D_256D(t1);
    x4 = GET_HI_128D_256D(t3);

    // t2 (ymm13) =
    // [
    //   P1010(S) + P1111(S)
    //   P1010(I) + P1111(I)
    //   P1010(J) + P1111(J)
    //   P1010(K) + P1111(K)
    // ]
    t2 = SET_256D_FROM_128D(x2,x4); 

    // C11 (ymm3) += t2 (ymm13)
    c11 = _mm256_add_pd(c11,t2);




    // COMPUTE OFF DIAGONAL CONTRIBUTIONS TO PRODUCT

    // Permute A
    //   a00 (ymm4) = [ A10(S) A11(S) A00(S) A01(S) ]
    //   a01 (ymm5) = [ A10(I) A11(I) A00(I) A01(I) ]
    //   a10 (ymm6) = [ A10(J) A11(J) A00(J) A01(J) ]
    //   a11 (ymm7) = [ A10(K) A11(K) A00(K) A01(K) ]

    // Use X1 (xmm12) and X2 (xmm13) to store the low and high bits

    x1 = GET_HI_128D_256D(a00); x2 = GET_LO_128D_256D(a00); 
      a00 = SET_256D_FROM_128D(x1,x2);
    x1 = GET_HI_128D_256D(a10); x2 = GET_LO_128D_256D(a10); 
      a10 = SET_256D_FROM_128D(x1,x2);
    x1 = GET_HI_128D_256D(a01); x2 = GET_LO_128D_256D(a01); 
      a01 = SET_256D_FROM_128D(x1,x2);
    x1 = GET_HI_128D_256D(a11); x2 = GET_LO_128D_256D(a11); 
      a11 = SET_256D_FROM_128D(x1,x2);

    // Compute product in T (ymm12-ymm15)
    // t1 (ymm12) = [ P1000(S) P1101(S) P0010(S) P0111(S) ]
    // t2 (ymm13) = [ P1000(I) P1101(I) P0010(I) P0111(I) ]
    // t3 (ymm14) = [ P1000(J) P1101(J) P0010(J) P0111(J) ]
    // t4 (ymm15) = [ P1000(K) P1101(K) P0010(K) P0111(K) ]
    // Multiplication requires no register stratch?
    MULD4Q_NN(a00,a01,a10,a11,b00,b01,b10,b11,t1,t2,t3,t4);


    // t1 (ymm12) =            t3 (ymm14) = 
    // [                       [
    //   P1000(S) + P1101(S)     P1000(J) + P1101(J)
    //   P1000(I) + P1101(I)     P1000(K) + P1101(K)
    //   P0010(S) + P0111(S)     P0010(J) + P0111(J)
    //   P0010(I) + P0111(I)     P0010(K) + P0111(K)
    // ]                       ]
    t1 = _mm256_hadd_pd(t1,t2);
    t3 = _mm256_hadd_pd(t3,t4);


    // XXX: Is there a better way to set the high bits of one
    // register to the low bits of another?


    // x2 (xmm13) =            x4 (xmm15) =
    // [                       [
    //   P1000(S) + P1101(S)     P1000(J) + P1101(J)
    //   P1000(I) + P1101(I)     P1000(K) + P1101(K)
    // ]                       ]
    x2 = GET_LO_128D_256D(t1);
    x4 = GET_LO_128D_256D(t3);


    // t2 (ymm13) =
    // [
    //   P1000(S) + P1101(S)
    //   P1000(I) + P1101(I)
    //   P1000(J) + P1101(J)
    //   P1000(K) + P1101(K)
    // ]
    t2 = SET_256D_FROM_128D(x2,x4); 
    
    // C10 (ymm1) += t2 (ymm13)
    c10 = _mm256_add_pd(c10,t2);


    // x2 (xmm13) =           x4 (xmm15) =
    // [                      [
    //   P0010(S) + P0111(S)    P0010(J) + P0111(J)
    //   P0010(I) + P0111(I)    P0010(K) + P0111(K)
    // ]                      ]
    x2 = GET_HI_128D_256D(t1);
    x4 = GET_HI_128D_256D(t3);

    // t2 (ymm13) =
    // [
    //   P0010(S) + P0111(S)
    //   P0010(I) + P0111(I)
    //   P0010(J) + P0111(J)
    //   P0010(K) + P0111(K)
    // ]
    t2 = SET_256D_FROM_128D(x2,x4); 

    // C01 (ymm2) += t2 (ymm13)
    c01 = _mm256_add_pd(c01,t2);
    

    k--;
  } while( k > 0 );


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
            Kern(BETA,nIII,nJJJ,nK,smallA,BL1,smallC,LDC);
     
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

