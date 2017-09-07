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
void HBLAS_GEMM(const char TRANSA, const char TRANSB, const HAXX_INT M, 
  const HAXX_INT N, const HAXX_INT K, const _AlphaF ALPHA, _AMatF * const A, 
  const HAXX_INT LDA, _BMatF * const B, const HAXX_INT LDB, 
  const _BetaF BETA, quaternion<_F> * const C, const HAXX_INT LDC){


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

  quaternion<_F> *CCol; 
  _BMatF *BCol, *BRow;
  _AMatF *ACol;

  // If ALPHA is zero
  if( AlphaIsZero ) {
    if( BetaIsZero )
      for( j = 0; j < N; ++j ) {
        CCol = C + j*LDC;
        for( i = 0; i < M; ++i ) CCol[i] = 0.;
      }
    else
      for( j = 0; j < N; ++j ) {
        CCol = C + j*LDC;
        for( i = 0; i < M; ++i ) CCol[i] = BETA * CCol[i];
      }

    return; // Nothing more to to
  }

  if( NOTB ) {
    if( NOTA ) {
      for( j = 0; j < N; ++j ) {
        CCol = C + j*LDC;
        BCol = B + j*LDB;

        if( BetaIsZero )
          for( i = 0; i < M; ++i ) CCol[i] = 0.;
        else if( not BetaIsOne )
          for( i = 0; i < M; ++i ) CCol[i] = BETA * CCol[i]; 

        // FIXME: Need to figure out a way to cache something in the
        // outer loop
        for( l = 0; l < K; ++l ) {
          ACol = A + l*LDA;
          for( i = 0; i < M; ++i ) CCol[i] += ALPHA * ACol[i] * BCol[l];
        }
      } // end loop-j
    } else if ( CONJA ) { // end NOTA
      for( j = 0; j < N; ++j ) {
        CCol = C + j*LDC;
        BCol = B + j*LDB;

        for( i = 0; i < M; ++i ) { 
          ACol = A + i*LDA;
          htemp = 0.;
          for( l = 0; l < K; ++l ) htemp += SmartConj(ACol[l]) * BCol[l];
  
          if( BetaIsZero ) CCol[i] = ALPHA*htemp;
          else             CCol[i] = ALPHA*htemp + BETA*CCol[i];
        } // end loop-i
      } // end loop-j
    } else { // end CONJA
      // Implies TRANSA == 'T'
      for( j = 0; j < N; ++j ) {
        CCol = C + j*LDC;
        BCol = B + j*LDB;

        for( i = 0; i < M; ++i ) { 
          ACol = A + i*LDA;
          htemp = 0.;
          for( l = 0; l < K; ++l ) htemp += ACol[l] * BCol[l];
  
          if( BetaIsZero ) CCol[i] = ALPHA*htemp;
          else             CCol[i] = ALPHA*htemp + BETA*CCol[i];
        } // end loop-i
      } // end loop-j
    } // end TRANSA 
  } else if( NOTA ) { // end NOTB
    if( CONJB ) {
      for( j = 0; j < N; ++j ) {
        CCol = C + j*LDC;

        if( BetaIsZero )
          for( i = 0; i < M; ++i ) CCol[i] = 0.;
        else if( not BetaIsOne )
          for( i = 0; i < M; ++i ) CCol[i] = BETA * CCol[i];

        // FIXME: Need to figure out a way to cache something in the
        // outer loop
        for( l = 0; l < K; ++l ) {
          ACol = A + l*LDA;
          BCol = B + l*LDB;
          for( i = 0; i < M; ++i ) CCol[i] += ALPHA * ACol[i] * SmartConj(BCol[j]);
        }
      } // end loop-j
    } else { // end CONJB
      // Implies TRANSB == 'T'
      for( j = 0; j < N; ++j ) {
        CCol = C + j*LDC;

        if( BetaIsZero )
          for( i = 0; i < M; ++i ) CCol[i] = 0.;
        else if( not BetaIsOne )
          for( i = 0; i < M; ++i ) CCol[i] = BETA * CCol[i];

        // FIXME: Need to figure out a way to cache something in the
        // outer loop
        for( l = 0; l < K; ++l ) {
          ACol = A + l*LDA;
          BCol = B + l*LDB;
          for( i = 0; i < M; ++i ) CCol[i] += ALPHA * ACol[i] * BCol[j];
        }
      } // end loop-j
    } // end TRANSB
  } else if( CONJA ) { // end NOTA
    if( CONJB ) {
      for( j = 0; j < N; ++j ) { 
        CCol = C + j*LDC;
        BRow = B + j;
        for( i = 0; i < M; ++i ) { 
          ACol = A + i*LDA;
          htemp = 0.;
          for( l = 0; l < K; ++l ) htemp += SmartConj(ACol[l]) * SmartConj(BRow[l*LDB]);

          if( BetaIsZero ) CCol[i] = ALPHA*htemp;
          else             CCol[i] = ALPHA*htemp + BETA*CCol[i];
        } // end loop-i
      } // end loop-j
    } else { // end CONJB
      for( j = 0; j < N; ++j ) { 
        CCol = C + j*LDC;
        BRow = B + j;
        for( i = 0; i < M; ++i ) { 
          ACol = A + i*LDA;
          htemp = 0.;
          for( l = 0; l < K; ++l ) htemp += SmartConj(ACol[l]) * BRow[l*LDB];

          if( BetaIsZero ) CCol[i] = ALPHA*htemp;
          else             CCol[i] = ALPHA*htemp + BETA*CCol[i];
        } // end loop-i
      } // end loop-j
    } // end TRANSB
  } else { // end CONJA
    // Implies TRANSA == 'T'
    if( CONJB ) {
      for( j = 0; j < N; ++j ) { 
        CCol = C + j*LDC;
        BRow = B + j;
        for( i = 0; i < M; ++i ) { 
          ACol = A + i*LDA;
          htemp = 0.;
          for( l = 0; l < K; ++l ) htemp += ACol[l] * SmartConj(BRow[l*LDB]);

          if( BetaIsZero ) CCol[i] = ALPHA*htemp;
          else             CCol[i] = ALPHA*htemp + BETA*CCol[i];
        } // end loop-i
      } // end loop-j
    } else { // end CONJB
      for( j = 0; j < N; ++j ) { 
        CCol = C + j*LDC;
        BRow = B + j;
        for( i = 0; i < M; ++i ) { 
          ACol = A + i*LDA;
          htemp = 0.;
          for( l = 0; l < K; ++l ) htemp += ACol[l] * BRow[l*LDB];

          if( BetaIsZero ) CCol[i] = ALPHA*htemp;
          else             CCol[i] = ALPHA*htemp + BETA*CCol[i];
        } // end loop-i
      } // end loop-j
    } // end TRANSB
  }
};

}; // namespace HAXX

#endif

