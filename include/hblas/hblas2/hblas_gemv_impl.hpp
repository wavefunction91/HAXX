/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_GEMV_IMPL_HPP
#define __INCLUDED_HBLAS_GEMV_IMPL_HPP

#include "hblas/hblas2_def.hpp"

namespace HAXX {

/**
 *  Written by DBWY (4/2017)
 *
 *  Based on the BLAS implementation of ZGEMV by
 *    Jack Dongarra   (Argonne)
 *    Jeremy Du Croz  (NAG)
 *    Sven Hammarling (NAG)
 *    Richard Hanson  (Sandia)
 *
 *  Performs one of the matrix - vector operations
 *
 *  \f$ y \in \mathbb{H} \qquad x,\alpha,A \in \mathbb{R},\mathbb{C},\mathbb{H}\f$
 *
 *  TRANS = 'N'
 *
 *  \f$ y_i = \alpha A_{ij} x_j + \beta y_i\f$
 *
 *  TRANS = 'T'
 *
 *  \f$ y_i = \alpha A_{ji} x_j + \beta y_i\f$
 *
 *  TRANS = 'C'
 *
 *  \f$ y_i = \alpha A^*_{ji} x_j + \beta y_i\f$
 *
 *
 */
// FIXME: In this implementaion, it has been implied that scalars
// will always multiply from the left. Should generalize in such a
// was to allow flexibility in both ALPHA and BETA
template <typename _F,typename _MatF, typename _VecF, typename _AlphaF, 
  typename _BetaF>
void HBLAS_GEMV(char TRANS, HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, 
  _MatF *A, HAXX_INT LDA, _VecF *X, HAXX_INT INCX, _BetaF BETA, 
  quaternion<_F> *Y, HAXX_INT INCY) {

  // FIXME: Should write functions to check quaternion equality
  //   to real / complex
  if( M == 0 or N == 0 or (ALPHA == _AlphaF(0.) and BETA == _BetaF(0.)) )
    return;

  // FIXME: The original BLAS implementaion has logic to handle
  //   negative strides. See further comments.
  assert( INCX > 0 );
  assert( INCY > 0 );

  bool NOCONJ = TRANS == 'T';

  HAXX_INT LENX = N, LENY = M;
  if( TRANS != 'N' ) {
    LENX = M;
    LENY = N;
  }

  // Indexing variables
  HAXX_INT i, j, iy, jx, ix, jy;

  // FIXME: These parameters are effected in the orignal BLAS
  //   implementaion by negative strides
  HAXX_INT KX = 0, KY = 0;
  quaternion<_F> htemp1, htemp2;

  // Form y = beta * y
  if( BETA != _BetaF(1.) ) {

    if( INCY == 1 ) {

      if( BETA == _BetaF(0.) ) 
        // XXX: STL call to fill_n?
        for( i = 0; i < LENY; ++i ) Y[i] = 0.; 
      else 
        // XXX: Why not use a level 1 BLAS call here?
        for( i = 0; i < LENY; ++i ) Y[i] = BETA * Y[i];

    } else { // end INCY == 1

      iy = KY;
      if( BETA == _BetaF(0.) ) 
        // XXX: STL call to fill_n
        for( i = 0; i < LENY; ++i, iy += INCY ) Y[iy] = 0.; 
      else
        // XXX: Why not use a level 1 BLAS call here?
        for( i = 0; i < LENY; ++i, iy += INCY ) Y[iy] = BETA * Y[iy]; 
      
    } // end INCY != 1

  } // end BETA == 1


  // Return if ALPHA == 0 as there's nothing more to  do
  if( ALPHA == _AlphaF(0.) ) return;

  if( TRANS == 'N' ) {

    // Form y_i = alpha * A_{ij} * x_j + beta * y_i

    jx = KX;
    if( INCY == 1 ) {

      for ( j = 0; j < N; ++j, jx += INCX ) {

        // alpha * A_{ij} * x_j = 
        //   [alpha,A_{ij}] * x_j + A_{ij} * (alpha * x_j)
        htemp1 = ALPHA * X[jx];
        for( i = 0; i < M; ++i ) {

          // y_j = y_j + A_{ij} * (alpha * x_j)
          Y[i] += A[RANK2_INDX(i,j,LDA)] * htemp1;

          // FIXME: This will probably kill vectorization
          //   Maybe just deal with the explicit left multiplication
          //   of Alpha?
          // y_j = y_j + [alpha,A_{ij}] * x_j
          Y[i] += comm(ALPHA,A[RANK2_INDX(i,j,LDA)]) * X[jx];

        } // end i loop

      } // end j loop

    } else { // end INCY == 1

      for ( j = 0; j < N; ++j, jx += INCX ) {

        // alpha * A_{ij} * x_j = 
        //   [alpha,A_{ij}] * x_j + A_{ij} * (alpha * x_j)
        htemp1 = ALPHA * X[jx];
        for( i = 0, iy = KY; i < M; ++i, iy += INCY ) {

          // y_j = y_j + A_{ij} * (alpha * x_j)
          Y[iy] += A[RANK2_INDX(i,j,LDA)] * htemp1;

          // FIXME: This will probably kill vectorization
          //   Maybe just deal with the explicit left multiplication
          //   of Alpha?
          // y_j = y_j + [alpha,A_{ij}] * x_j
          Y[iy] += comm(ALPHA,A[RANK2_INDX(i,j,LDA)]) * X[jx];

        } // end i loop

      } // end j loop

    } // end INCY != 1

  } else { // end TRANS == 'N'

    // form y_i = alpha * A_{ji} x_j + beta * y_i 
    //
    //   or 
    //
    // form y_i = alpha * A^*_{ji} x_j + beta * y_i 

    jy = KY;
    if( INCX == 1 ) {

      for( j = 0; j < N; ++j, jy += INCY ) {

        // XXX: Why not use a level 1 BLAS call here?
        htemp1 = 0.;
        if( NOCONJ )
          for( i = 0; i < M; ++i) 
            htemp1 += A[RANK2_INDX(i,j,LDA)] * X[i];
        else
          for( i = 0; i < M; ++i) 
            htemp1 += conj(A[RANK2_INDX(i,j,LDA)]) * X[i];

        Y[jy] += ALPHA * htemp1;

      } // end j loop

    } else { // end INCX == 1

      for( j = 0; j < N; ++j, jy += INCY ) {

        // XXX: Why not use a level 1 BLAS call here?
        htemp1 = 0.;
        if( NOCONJ )
          for( i = 0, ix = KX; i < M; ++i, ix += INCX) 
            htemp1 += A[RANK2_INDX(i,j,LDA)] * X[ix];
        else
          for( i = 0, ix = KX; i < M; ++i, ix += INCX) 
            htemp1 += conj(A[RANK2_INDX(i,j,LDA)]) * X[ix];

        Y[jy] += ALPHA * htemp1;

      } // end j loop

    } // end INCX != 1

  } // end TRANS != 'N'

}; // end GEMV




// Specialization for REAL Alpha
template <typename _F,typename _MatF, typename _VecF, typename _BetaF>
void HBLAS_GEMV(char TRANS, HAXX_INT M, HAXX_INT N, _F ALPHA, 
  _MatF *A, HAXX_INT LDA, _VecF *X, HAXX_INT INCX, _BetaF BETA, 
  quaternion<_F> *Y, HAXX_INT INCY) {

  // FIXME: Should write functions to check quaternion equality
  //   to real / complex
  if( M == 0 or N == 0 or (ALPHA == 0. and BETA == _BetaF(0.)) )
    return;

  // FIXME: The original BLAS implementaion has logic to handle
  //   negative strides. See further comments.
  assert( INCX > 0 );
  assert( INCY > 0 );

  bool NOCONJ = TRANS == 'T';

  HAXX_INT LENX = N, LENY = M;
  if( TRANS != 'N' ) {
    LENX = M;
    LENY = N;
  }

  // Indexing variables
  HAXX_INT i, j, iy, jx, ix, jy;

  // FIXME: These parameters are effected in the orignal BLAS
  //   implementaion by negative strides
  HAXX_INT KX = 0, KY = 0;
  quaternion<_F> htemp1, htemp2;

  // Form y = beta * y
  if( BETA != _BetaF(1.) ) {

    if( INCY == 1 ) {

      if( BETA == _BetaF(0.) ) 
        // XXX: STL call to fill_n?
        for( i = 0; i < LENY; ++i ) Y[i] = 0.; 
      else 
        // XXX: Why not use a level 1 BLAS call here?
        for( i = 0; i < LENY; ++i ) Y[i] = BETA * Y[i];

    } else { // end INCY == 1

      iy = KY;
      if( BETA == _BetaF(0.) ) 
        // XXX: STL call to fill_n
        for( i = 0; i < LENY; ++i, iy += INCY ) Y[iy] = 0.; 
      else
        // XXX: Why not use a level 1 BLAS call here?
        for( i = 0; i < LENY; ++i, iy += INCY ) Y[iy] = BETA * Y[iy]; 
      
    } // end INCY != 1

  } // end BETA == 1


  // Return if ALPHA == 0 as there's nothing more to  do
  if( ALPHA == 0. ) return;

  if( TRANS == 'N' ) {

    // Form y_i = alpha * A_{ij} * x_j + beta * y_i

    jx = KX;
    if( INCY == 1 ) {

      for ( j = 0; j < N; ++j, jx += INCX ) {

        // y_j = y_j + A_{ij} * (alpha * x_j)
        htemp1 = ALPHA * X[jx];
        for( i = 0; i < M; ++i ) 
          Y[i] += A[RANK2_INDX(i,j,LDA)] * htemp1;

      } // end j loop

    } else { // end INCY == 1

      for ( j = 0; j < N; ++j, jx += INCX ) {

        // y_j = y_j + A_{ij} * (alpha * x_j)
        htemp1 = ALPHA * X[jx];
        for( i = 0, iy = KY; i < M; ++i, iy += INCY ) 
          Y[iy] += A[RANK2_INDX(i,j,LDA)] * htemp1;


      } // end j loop

    } // end INCY != 1

  } else { // end TRANS == 'N'

    // form y_i = alpha * A_{ji} x_j + beta * y_i 
    //
    //   or 
    //
    // form y_i = alpha * A^*_{ji} x_j + beta * y_i 

    jy = KY;
    if( INCX == 1 ) {

      for( j = 0; j < N; ++j, jy += INCY ) {

        // XXX: Why not use a level 1 BLAS call here?
        htemp1 = 0.;
        if( NOCONJ )
          for( i = 0; i < M; ++i) 
            htemp1 += A[RANK2_INDX(i,j,LDA)] * X[i];
        else
          for( i = 0; i < M; ++i) 
            htemp1 += conj(A[RANK2_INDX(i,j,LDA)]) * X[i];

        Y[jy] += ALPHA * htemp1;

      } // end j loop

    } else { // end INCX == 1

      for( j = 0; j < N; ++j, jy += INCY ) {

        // XXX: Why not use a level 1 BLAS call here?
        htemp1 = 0.;
        if( NOCONJ )
          for( i = 0, ix = KX; i < M; ++i, ix += INCX) 
            htemp1 += A[RANK2_INDX(i,j,LDA)] * X[ix];
        else
          for( i = 0, ix = KX; i < M; ++i, ix += INCX) 
            htemp1 += conj(A[RANK2_INDX(i,j,LDA)]) * X[ix];

        Y[jy] += ALPHA * htemp1;

      } // end j loop

    } // end INCX != 1

  } // end TRANS != 'N'

}; // end GEMV (REAL ALPHA Specialization)

}; // namespace HAXX

#endif
