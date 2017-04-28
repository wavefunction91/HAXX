/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_GEMM_IMPL_HPP
#define __INCLUDED_HBLAS_GEMM_IMPL_HPP

#include "hblas/hblas3_def.hpp"

namespace HAXX {

template <typename _F, typename _AMatF, typename _BMatF, typename _AlphaF, 
  typename _BetaF>
void HBLAS_GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT N, HAXX_INT K,
  _AlphaF ALPHA, _AMatF *A, HAXX_INT LDA, _BMatF *B, HAXX_INT LDB, 
  _BetaF BETA, quaternion<_F> *C, HAXX_INT LDC){

  bool NOTA  = TRANSA == 'N';
  bool NOTB  = TRANSB == 'N';
  bool CONJA = TRANSA == 'C';
  bool CONJB = TRANSB == 'C';

  HAXX_INT NROWA, NCOLA, NROWB;

  if( NOTA ) {
    NROWA = M;
    NCOLA = K;
  } else {
    NROWA = K;
    NCOLA = M;
  }

  if( NOTB ) NROWB = K;
  else       NROWB = N;

  // FIXME: Add asserts for input sanity check

  _AlphaF AlphaZero = _AlphaF(0.);
  _AlphaF AlphaOne  = _AlphaF(1.);

  _BetaF BetaZero = _BetaF(0.);
  _BetaF BetaOne  = _BetaF(1.);


  bool BetaIsZero = BETA == BetaZero;
  bool BetaIsOne  = BETA == BetaOne;
  bool AlphaIsZero = ALPHA == AlphaZero;
  bool AlphaIsOne  = ALPHA == AlphaOne;

  // Quick return
  if( (M == 0) or (N == 0) or
    ( (AlphaIsZero or (K == 0)) and BetaIsOne )) return;

  HAXX_INT i,j,l;
  quaternion<_F> htemp;

  // If ALPHA is zero
  if( AlphaIsZero ) {
    if( BetaIsZero )
      for( j = 0; j < N; ++j)
      for( i = 0; i < M; ++i)
        C[RANK2_INDX(i,j,LDC)] = 0.;
    else
      for( j = 0; j < N; ++j)
      for( i = 0; i < M; ++i) 
        C[RANK2_INDX(i,j,LDC)] = BETA * C[RANK2_INDX(i,j,LDC)];

    return; // Nothing more to to
  }

  if( NOTB ) {
    if( NOTA ) {
      for( j = 0; j < N; ++j ) {
        if( BetaIsZero )
          for( i = 0; i < M; ++i) C[RANK2_INDX(i,j,LDC)] = 0.;
        else if( not BetaIsOne )
          for( i = 0; i < M; ++i) 
            C[RANK2_INDX(i,j,LDC)] = BETA * C[RANK2_INDX(i,j,LDC)];

        // FIXME: Need to figure out a way to cache something in the
        // outer loop
        for( l = 0; l < K; ++l )
        for( i = 0; i < M; ++i )
          C[RANK2_INDX(i,j,LDC)] += 
            ALPHA * A[RANK2_INDX(i,l,LDA)] * B[RANK2_INDX(l,j,LDB)];
      } // end loop-j
    } else if ( CONJA ) { // end NOTA
      for( j = 0; j < N; ++j ) 
      for( i = 0; i < M; ++i ) { 
        htemp = 0.;
        for( l = 0; l < K; ++l ) 
          htemp += conj(A[RANK2_INDX(l,i,LDA)]) * B[RANK2_INDX(l,j,LDB)];

        if( BetaIsZero ) 
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp;
        else
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp + BETA*C[RANK2_INDX(i,j,LDC)];

      } // end loop-ij
    } else { // end CONJA
      // Implies TRANSA == 'T'
      for( j = 0; j < N; ++j ) 
      for( i = 0; i < M; ++i ) { 
        htemp = 0.;
        for( l = 0; l < K; ++l ) 
          htemp += A[RANK2_INDX(l,i,LDA)] * B[RANK2_INDX(l,j,LDB)];

        if( BetaIsZero ) 
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp;
        else
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp + BETA*C[RANK2_INDX(i,j,LDC)];

      } // end loop-ij
    } // end TRANSA 
  } else if( NOTA ) { // end NOTB
    if( CONJB ) {
      for( j = 0; j < N; ++j ) {
        if( BetaIsZero )
          for( i = 0; i < M; ++i) C[RANK2_INDX(i,j,LDC)] = 0.;
        else if( not BetaIsOne )
          for( i = 0; i < M; ++i) 
            C[RANK2_INDX(i,j,LDC)] = BETA * C[RANK2_INDX(i,j,LDC)];

        // FIXME: Need to figure out a way to cache something in the
        // outer loop
        for( l = 0; l < K; ++l )
        for( i = 0; i < M; ++i )
          C[RANK2_INDX(i,j,LDC)] += 
            ALPHA * A[RANK2_INDX(i,l,LDA)] * conj(B[RANK2_INDX(j,l,LDB)]);
      } // end loop-j
    } else { // end CONJB
      // Implies TRANSB == 'T'
      for( j = 0; j < N; ++j ) {
        if( BetaIsZero )
          for( i = 0; i < M; ++i) C[RANK2_INDX(i,j,LDC)] = 0.;
        else if( not BetaIsOne )
          for( i = 0; i < M; ++i) 
            C[RANK2_INDX(i,j,LDC)] = BETA * C[RANK2_INDX(i,j,LDC)];

        // FIXME: Need to figure out a way to cache something in the
        // outer loop
        for( l = 0; l < K; ++l )
        for( i = 0; i < M; ++i )
          C[RANK2_INDX(i,j,LDC)] += 
            ALPHA * A[RANK2_INDX(i,l,LDA)] * B[RANK2_INDX(j,l,LDB)];
      } // end loop-j
    } // end TRANSB
  } else if( CONJA ) { // end NOTA
    if( CONJB ) {
      for( j = 0; j < N; ++j ) 
      for( i = 0; i < M; ++i ) { 
        htemp = 0.;
        for( l = 0; l < K; ++l ) 
          htemp += conj(A[RANK2_INDX(l,i,LDA)]) * conj(B[RANK2_INDX(j,l,LDB)]);

        if( BetaIsZero ) 
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp;
        else
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp + BETA*C[RANK2_INDX(i,j,LDC)];
      } // end loop-ij
    } else { // end CONJB
      for( j = 0; j < N; ++j ) 
      for( i = 0; i < M; ++i ) { 
        htemp = 0.;
        for( l = 0; l < K; ++l ) 
          htemp += conj(A[RANK2_INDX(l,i,LDA)]) * B[RANK2_INDX(j,l,LDB)];

        if( BetaIsZero ) 
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp;
        else
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp + BETA*C[RANK2_INDX(i,j,LDC)];
      } // end loop-ij
    } // end TRANSB
  } else { // end CONJA
    // Implies TRANSA == 'T'
    if( CONJB ) {
      for( j = 0; j < N; ++j ) 
      for( i = 0; i < M; ++i ) { 
        htemp = 0.;
        for( l = 0; l < K; ++l ) 
          htemp += A[RANK2_INDX(l,i,LDA)] * conj(B[RANK2_INDX(j,l,LDB)]);

        if( BetaIsZero ) 
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp;
        else
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp + BETA*C[RANK2_INDX(i,j,LDC)];
      } // end loop-ij
    } else { // end CONJB
      for( j = 0; j < N; ++j ) 
      for( i = 0; i < M; ++i ) { 
        htemp = 0.;
        for( l = 0; l < K; ++l ) 
          htemp += A[RANK2_INDX(l,i,LDA)] * B[RANK2_INDX(j,l,LDB)];

        if( BetaIsZero ) 
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp;
        else
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp + BETA*C[RANK2_INDX(i,j,LDC)];
      } // end loop-ij
    } // end TRANSB
  }
};



template <typename _F, typename _AMatF, typename _BMatF, typename _BetaF>
void HBLAS_GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT N, HAXX_INT K,
  _F ALPHA, _AMatF *A, HAXX_INT LDA, _BMatF *B, HAXX_INT LDB, 
  _BetaF BETA, quaternion<_F> *C, HAXX_INT LDC){

  bool NOTA  = TRANSA == 'N';
  bool NOTB  = TRANSB == 'N';
  bool CONJA = TRANSA == 'C';
  bool CONJB = TRANSB == 'C';

  HAXX_INT NROWA, NCOLA, NROWB;

  if( NOTA ) {
    NROWA = M;
    NCOLA = K;
  } else {
    NROWA = K;
    NCOLA = M;
  }

  if( NOTB ) NROWB = K;
  else       NROWB = N;

  // FIXME: Add asserts for input sanity check

  _F AlphaZero = _F(0.);
  _F AlphaOne  = _F(1.);

  _BetaF BetaZero = _BetaF(0.);
  _BetaF BetaOne  = _BetaF(1.);


  bool BetaIsZero = BETA == BetaZero;
  bool BetaIsOne  = BETA == BetaOne;
  bool AlphaIsZero = ALPHA == AlphaZero;
  bool AlphaIsOne  = ALPHA == AlphaOne;

  // Quick return
  if( (M == 0) or (N == 0) or
    ( (AlphaIsZero or (K == 0)) and BetaIsOne )) return;

  HAXX_INT i,j,l;
  quaternion<_F> htemp;

  // If ALPHA is zero
  if( AlphaIsZero ) {
    if( BetaIsZero )
      for( j = 0; j < N; ++j)
      for( i = 0; i < M; ++i)
        C[RANK2_INDX(i,j,LDC)] = 0.;
    else
      for( j = 0; j < N; ++j)
      for( i = 0; i < M; ++i) 
        C[RANK2_INDX(i,j,LDC)] = BETA * C[RANK2_INDX(i,j,LDC)];

    return; // Nothing more to to
  }

  if( NOTB ) {
    if( NOTA ) {
      for( j = 0; j < N; ++j ) {
        if( BetaIsZero )
          for( i = 0; i < M; ++i) C[RANK2_INDX(i,j,LDC)] = 0.;
        else if( not BetaIsOne )
          for( i = 0; i < M; ++i) 
            C[RANK2_INDX(i,j,LDC)] = BETA * C[RANK2_INDX(i,j,LDC)];

        for( l = 0; l < K; ++l ) {
          htemp = ALPHA * B[RANK2_INDX(l,j,LDB)];
          for( i = 0; i < M; ++i )
            C[RANK2_INDX(i,j,LDC)] += A[RANK2_INDX(i,l,LDA)] * htemp;
        } // end loop-l
      } // end loop-j
    } else if ( CONJA ) { // end NOTA
      for( j = 0; j < N; ++j ) 
      for( i = 0; i < M; ++i ) { 
        htemp = 0.;
        for( l = 0; l < K; ++l ) 
          htemp += conj(A[RANK2_INDX(l,i,LDA)]) * B[RANK2_INDX(l,j,LDB)];

        if( BetaIsZero ) 
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp;
        else
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp + BETA*C[RANK2_INDX(i,j,LDC)];

      } // end loop-ij
    } else { // end CONJA
      // Implies TRANSA == 'T'
      for( j = 0; j < N; ++j ) 
      for( i = 0; i < M; ++i ) { 
        htemp = 0.;
        for( l = 0; l < K; ++l ) 
          htemp += A[RANK2_INDX(l,i,LDA)] * B[RANK2_INDX(l,j,LDB)];

        if( BetaIsZero ) 
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp;
        else
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp + BETA*C[RANK2_INDX(i,j,LDC)];

      } // end loop-ij
    } // end TRANSA 
  } else if( NOTA ) { // end NOTB
    if( CONJB ) {
      for( j = 0; j < N; ++j ) {
        if( BetaIsZero )
          for( i = 0; i < M; ++i) C[RANK2_INDX(i,j,LDC)] = 0.;
        else if( not BetaIsOne )
          for( i = 0; i < M; ++i) 
            C[RANK2_INDX(i,j,LDC)] = BETA * C[RANK2_INDX(i,j,LDC)];

        for( l = 0; l < K; ++l ) {
          htemp = ALPHA * conj(B[RANK2_INDX(j,l,LDB)]);
          for( i = 0; i < M; ++i )
            C[RANK2_INDX(i,j,LDC)] += A[RANK2_INDX(i,l,LDA)] * htemp;
        } // end loop-l
      } // end loop-j
    } else { // end CONJB
      // Implies TRANSB == 'T'
      for( j = 0; j < N; ++j ) {
        if( BetaIsZero )
          for( i = 0; i < M; ++i) C[RANK2_INDX(i,j,LDC)] = 0.;
        else if( not BetaIsOne )
          for( i = 0; i < M; ++i) 
            C[RANK2_INDX(i,j,LDC)] = BETA * C[RANK2_INDX(i,j,LDC)];

        for( l = 0; l < K; ++l ) {
          htemp = ALPHA * B[RANK2_INDX(j,l,LDB)];
          for( i = 0; i < M; ++i )
            C[RANK2_INDX(i,j,LDC)] += A[RANK2_INDX(i,l,LDA)] * htemp;
        } // end loop-l
      } // end loop-j
    } // end TRANSB
  } else if( CONJA ) { // end NOTA
    if( CONJB ) {
      for( j = 0; j < N; ++j ) 
      for( i = 0; i < M; ++i ) { 
        htemp = 0.;
        for( l = 0; l < K; ++l ) 
          htemp += conj(A[RANK2_INDX(l,i,LDA)]) * conj(B[RANK2_INDX(j,l,LDB)]);

        if( BetaIsZero ) 
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp;
        else
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp + BETA*C[RANK2_INDX(i,j,LDC)];
      } // end loop-ij
    } else { // end CONJB
      for( j = 0; j < N; ++j ) 
      for( i = 0; i < M; ++i ) { 
        htemp = 0.;
        for( l = 0; l < K; ++l ) 
          htemp += conj(A[RANK2_INDX(l,i,LDA)]) * B[RANK2_INDX(j,l,LDB)];

        if( BetaIsZero ) 
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp;
        else
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp + BETA*C[RANK2_INDX(i,j,LDC)];
      } // end loop-ij
    } // end TRANSB
  } else { // end CONJA
    // Implies TRANSA == 'T'
    if( CONJB ) {
      for( j = 0; j < N; ++j ) 
      for( i = 0; i < M; ++i ) { 
        htemp = 0.;
        for( l = 0; l < K; ++l ) 
          htemp += A[RANK2_INDX(l,i,LDA)] * conj(B[RANK2_INDX(j,l,LDB)]);

        if( BetaIsZero ) 
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp;
        else
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp + BETA*C[RANK2_INDX(i,j,LDC)];
      } // end loop-ij
    } else { // end CONJB
      for( j = 0; j < N; ++j ) 
      for( i = 0; i < M; ++i ) { 
        htemp = 0.;
        for( l = 0; l < K; ++l ) 
          htemp += A[RANK2_INDX(l,i,LDA)] * B[RANK2_INDX(j,l,LDB)];

        if( BetaIsZero ) 
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp;
        else
          C[RANK2_INDX(i,j,LDC)] = ALPHA*htemp + BETA*C[RANK2_INDX(i,j,LDC)];
      } // end loop-ij
    } // end TRANSB
  }
};

template<>
void HBLAS_GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT N, HAXX_INT K,
  double ALPHA, quaternion<double> *A, HAXX_INT LDA, quaternion<double> *B, HAXX_INT LDB, 
  double BETA, quaternion<double> *C, HAXX_INT LDC){

  std::cout << " IN SDD" << std::endl;
  hgemmdd_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,reinterpret_cast<double*>(A),
    &LDA,reinterpret_cast<double*>(B),&LDB,&BETA,
    reinterpret_cast<double*>(C),&LDC);

};

template<>
void HBLAS_GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT N, HAXX_INT K,
  double ALPHA, quaternion<double> *A, HAXX_INT LDA, quaternion<double> *B, HAXX_INT LDB, 
  std::complex<double> BETA, quaternion<double> *C, HAXX_INT LDC){

  std::cout << " IN SDZ" << std::endl;
  hgemmdz_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,reinterpret_cast<double*>(A),
    &LDA,reinterpret_cast<double*>(B),&LDB,reinterpret_cast<double*>(&BETA),
    reinterpret_cast<double*>(C),&LDC);

};

template<>
void HBLAS_GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT N, HAXX_INT K,
  double ALPHA, quaternion<double> *A, HAXX_INT LDA, quaternion<double> *B, HAXX_INT LDB, 
  quaternion<double> BETA, quaternion<double> *C, HAXX_INT LDC){

  std::cout << " IN SDH" << std::endl;
  hgemmdh_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,reinterpret_cast<double*>(A),
    &LDA,reinterpret_cast<double*>(B),&LDB,reinterpret_cast<double*>(&BETA),
    reinterpret_cast<double*>(C),&LDC);

};

template<>
void HBLAS_GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT N, HAXX_INT K,
  std::complex<double> ALPHA, quaternion<double> *A, HAXX_INT LDA, quaternion<double> *B, HAXX_INT LDB, 
  double BETA, quaternion<double> *C, HAXX_INT LDC){

  hgemmzd_(&TRANSA,&TRANSB,&M,&N,&K,reinterpret_cast<double*>(&ALPHA),
    reinterpret_cast<double*>(A),&LDA,reinterpret_cast<double*>(B),
    &LDB,&BETA, reinterpret_cast<double*>(C),&LDC);

};

template<>
void HBLAS_GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT N, HAXX_INT K,
  std::complex<double> ALPHA, quaternion<double> *A, HAXX_INT LDA, quaternion<double> *B, HAXX_INT LDB, 
  std::complex<double> BETA, quaternion<double> *C, HAXX_INT LDC){

  std::cout << " IN SZZ" << std::endl;
  hgemmzz_(&TRANSA,&TRANSB,&M,&N,&K,reinterpret_cast<double*>(&ALPHA),
    reinterpret_cast<double*>(A),&LDA,reinterpret_cast<double*>(B),
    &LDB,reinterpret_cast<double*>(&BETA), reinterpret_cast<double*>(C),&LDC);

};

template<>
void HBLAS_GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT N, HAXX_INT K,
  std::complex<double> ALPHA, quaternion<double> *A, HAXX_INT LDA, quaternion<double> *B, HAXX_INT LDB, 
  quaternion<double> BETA, quaternion<double> *C, HAXX_INT LDC){

  std::cout << " IN SZH" << std::endl;
  hgemmzh_(&TRANSA,&TRANSB,&M,&N,&K,reinterpret_cast<double*>(&ALPHA),
    reinterpret_cast<double*>(A),&LDA,reinterpret_cast<double*>(B),
    &LDB,reinterpret_cast<double*>(&BETA), reinterpret_cast<double*>(C),&LDC);

};

template<>
void HBLAS_GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT N, HAXX_INT K,
  quaternion<double> ALPHA, quaternion<double> *A, HAXX_INT LDA, quaternion<double> *B, HAXX_INT LDB, 
  double BETA, quaternion<double> *C, HAXX_INT LDC){

  hgemmhd_(&TRANSA,&TRANSB,&M,&N,&K,reinterpret_cast<double*>(&ALPHA),
    reinterpret_cast<double*>(A),&LDA,reinterpret_cast<double*>(B),
    &LDB,&BETA, reinterpret_cast<double*>(C),&LDC);

};

template<>
void HBLAS_GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT N, HAXX_INT K,
  quaternion<double> ALPHA, quaternion<double> *A, HAXX_INT LDA, quaternion<double> *B, HAXX_INT LDB, 
  std::complex<double> BETA, quaternion<double> *C, HAXX_INT LDC){

  hgemmhz_(&TRANSA,&TRANSB,&M,&N,&K,reinterpret_cast<double*>(&ALPHA),
    reinterpret_cast<double*>(A),&LDA,reinterpret_cast<double*>(B),
    &LDB,reinterpret_cast<double*>(&BETA), reinterpret_cast<double*>(C),&LDC);

};

template<>
void HBLAS_GEMM(char TRANSA, char TRANSB, HAXX_INT M, HAXX_INT N, HAXX_INT K,
  quaternion<double> ALPHA, quaternion<double> *A, HAXX_INT LDA, quaternion<double> *B, HAXX_INT LDB, 
  quaternion<double> BETA, quaternion<double> *C, HAXX_INT LDC){

  std::cout << " In SHH" << std::endl;
  hgemmhh_(&TRANSA,&TRANSB,&M,&N,&K,reinterpret_cast<double*>(&ALPHA),
    reinterpret_cast<double*>(A),&LDA,reinterpret_cast<double*>(B),
    &LDB,reinterpret_cast<double*>(&BETA), reinterpret_cast<double*>(C),&LDC);

};

}; // namespace HAXX

#endif

