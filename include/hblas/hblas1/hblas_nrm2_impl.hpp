/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_NRM2_IMPL_HPP
#define __INCLUDED_HBLAS_NRM2_IMPL_HPP

#include "hblas/hblas1_def.hpp"

namespace HAXX {

/**
 *  Written by DBWY (5/2017)
 *
 *  Based on the BLAS implementation of DZNRM2 by Sven Hammarling
 *    http://www.netlib.org/lapack/explore-html/d9/d19/dznrm2_8f.html 
 *
 *  \f$ r \in \mathbb{R}, \quad x \in \mathbb{H}^n \qquad r = x_i^* x_i
 */
template <typename _F>
_F HBLAS_NRM2( HAXX_INT N, quaternion<_F> *X, HAXX_INT INCX ) {

 if( N < 1 or INCX < 1 ) return _F(0.);
 else {
   _F scale = _F(0.);
   _F ssq   = _F(0.);
   _F temp;

   HAXX_INT IX = 0;

   for( IX = 0; IX < (N-1)*INCX + 1; IX += INCX ) {

     if( X[IX].real() != 0. ) {
       temp = std::abs(X[IX].real());
       if( scale < temp ) {
         ssq = 1. + ssq * (scale / temp) * (scale / temp);
         scale = temp;
       } else {
         ssq += (temp / scale) * (temp / scale);
       }
     }

     if( X[IX].imag_i() != 0. ) {
       temp = std::abs(X[IX].imag_i());
       if( scale < temp ) {
         ssq = 1. + ssq * (scale / temp) * (scale / temp);
         scale = temp;
       } else {
         ssq += (temp / scale) * (temp / scale);
       }
     }

     if( X[IX].imag_j() != 0. ) {
       temp = std::abs(X[IX].imag_j());
       if( scale < temp ) {
         ssq = 1. + ssq * (scale / temp) * (scale / temp);
         scale = temp;
       } else {
         ssq += (temp / scale) * (temp / scale);
       }
     }

     if( X[IX].imag_k() != 0. ) {
       temp = std::abs(X[IX].imag_k());
       if( scale < temp ) {
         ssq = 1. + ssq * (scale / temp) * (scale / temp);
         scale = temp;
       } else {
         ssq += (temp / scale) * (temp / scale);
       }
     }

   }
    
   return scale * std::sqrt(ssq);
 } 

};

/*
template<>
double HBLAS_NRM2( HAXX_INT N, quaternion<double> *X, HAXX_INT INCX ) {
  return dhnrm2_(&N,X,&INCX);
};
*/

}; // namespace HAXX

#endif
