/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_GERU_IMPL_HPP
#define __INCLUDED_HBLAS_GERU_IMPL_HPP

#include "hblas/hblas2_def.hpp"

namespace HAXX {

/**
 *  Written by DBWY (4/2017)
 *
 *  Based on the BLAS implementation of ZGERU by
 *    Jack Dongarra   (Argonne)
 *    Jeremy Du Croz  (NAG)
 *    Sven Hammarling (NAG)
 *    Richard Hanson  (Sandia)
 *
 *  Performs the rank 1 operation
 *
 *  \f$ A_{ij} = A_{ij} + \alpha x_i y_j \f$ 
 *
 *  \f$ A \in \mathbb{H} \qquad x,y,\alpha \in \mathbb{R},\mathbb{C},\mathbb{H}\f$
 */
// FIXME: In this implementaion, it has been implied that scalars
// will always multiply from the left. Should generalize in such a
// was to allow flexibility in ALPHA 
template <typename _F, typename _LeftVecF, typename _RightVecF, typename _AlphaF>
void HBLAS_GERU(HAXX_INT M, HAXX_INT N, _AlphaF ALPHA, _LeftVecF *X,
  HAXX_INT INCX, _RightVecF *Y, HAXX_INT INCY, quaternion<_F> *A, 
  HAXX_INT LDA){

  if( M == 0 or N == 0 or ALPHA == _AlphaF(0.)) return;

  // FIXME: The original BLAS implementaion has logic to handle
  //   negative strides. See further comments.
  assert( INCX > 0 );
  assert( INCY > 0 );


  HAXX_INT i, j, ix;
  
  // FIXME: This parameter is effected in the orignal BLAS
  //   implementaion by negative stride
  HAXX_INT JY = 0;

  quaternion<_F> htemp1;
 
  if( INCX == 1 ) {

    for( j = 0; j < N; ++j, JY += INCY )
      if( Y[JY] != _RightVecF(0.) ) {
//      htemp1 = ALPHA * Y[JY];
        for( i = 0; i < M; ++i ) {
          A[RANK2_INDX(i,j,LDA)] += ALPHA * X[i] * Y[JY];
        }
      }

  } else { // end INCX == 1

    // FIXME: This parameter is effected in the orignal BLAS
    //   implementaion by negative stride
    HAXX_INT KX = 0;

    for( j = 0; j < N; ++j, JY += INCY )
      if( Y[JY] != _RightVecF(0.) ) {
//      htemp1 = ALPHA * Y[JY];
        for( i = 0, ix = KX; i < M; ++i, ix += INCX ) {
          A[RANK2_INDX(i,j,LDA)] += ALPHA * X[ix] * Y[JY];
        }
      }
  } // end INCX != 1

}; // end GERU


template <>
void HBLAS_GERU(HAXX_INT M, HAXX_INT N, double ALPHA, quaternion<double> *X,
  HAXX_INT INCX, quaternion<double> *Y, HAXX_INT INCY, quaternion<double> *A, 
  HAXX_INT LDA){

  hgerud_(&M,&N,&ALPHA,reinterpret_cast<double*>(X),&INCX,
    reinterpret_cast<double*>(Y),&INCY,reinterpret_cast<double*>(A),
    &LDA);
}

template <>
void HBLAS_GERU(HAXX_INT M, HAXX_INT N, std::complex<double> ALPHA, quaternion<double> *X,
  HAXX_INT INCX, quaternion<double> *Y, HAXX_INT INCY, quaternion<double> *A, 
  HAXX_INT LDA){

  hgeruz_(&M,&N,reinterpret_cast<double*>(&ALPHA),reinterpret_cast<double*>(X),&INCX,
    reinterpret_cast<double*>(Y),&INCY,reinterpret_cast<double*>(A),
    &LDA);
}

template <>
void HBLAS_GERU(HAXX_INT M, HAXX_INT N, quaternion<double> ALPHA, quaternion<double> *X,
  HAXX_INT INCX, quaternion<double> *Y, HAXX_INT INCY, quaternion<double> *A, 
  HAXX_INT LDA){

  hgeruh_(&M,&N,reinterpret_cast<double*>(&ALPHA),reinterpret_cast<double*>(X),&INCX,
    reinterpret_cast<double*>(Y),&INCY,reinterpret_cast<double*>(A),
    &LDA);
}

}; // namespace HAXX

#endif
