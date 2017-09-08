/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#include "haxx.hpp"
#include "hblas/hblas1_impl.hpp"



// Preprocessor macros for quick CXX implementations

#define SWAP_CXX_IMPL(F) \
template \
void HBLAS_SWAP(const HAXX_INT, quaternion<F> * const , const HAXX_INT,\
  quaternion<F> * const , const HAXX_INT);

#define COPY_CXX_IMPL(F) \
template \
void HBLAS_COPY(const HAXX_INT, quaternion<F> * const , const HAXX_INT,\
  quaternion<F> * const , const HAXX_INT);



// Preprocessor macros for quick FORTRAN implementations

#define DOT_FORTRAN_IMPL(FUNCNAME,NAME,F)\
template<>\
quaternion<F> HBLAS_##FUNCNAME( HAXX_INT N, quaternion<F> *X, \
  HAXX_INT INCX, quaternion<F> *Y, HAXX_INT INCY){\
  \
  quaternion<F> htemp;\
  NAME##_(&htemp, &N, X, &INCX, Y, &INCY);\
  return htemp;\
  \
}

#define SCAL_FORTRAN_IMPL(NAME,F,ALPHAF)\
template<>\
void HBLAS_SCAL(const char SIDE, const HAXX_INT N, const ALPHAF ALPHA,\
  quaternion<F> * const X, const HAXX_INT INCX) {\
  \
  NAME##_(&SIDE,&N,&ALPHA,X,&INCX);\
  \
};


#define AXPY_FORTRAN_IMPL(NAME,F,XF,ALPHAF)\
template<> \
void HBLAS_AXPY(const char SIDE, const HAXX_INT N, const ALPHAF ALPHA,\
  XF * const X, const HAXX_INT INCX, quaternion<F> * const Y, \
  const HAXX_INT INCY) {\
  \
  NAME##_(&SIDE,&N,&ALPHA,X,&INCX,Y,&INCY);\
  \
};

namespace HAXX {

  SWAP_CXX_IMPL(double);
  
  SCAL_FORTRAN_IMPL(hscald,double,double);
  SCAL_FORTRAN_IMPL(hscalc,double,std::complex<double>);
  SCAL_FORTRAN_IMPL(hscalh,double,quaternion<double>);

  COPY_CXX_IMPL(double);

  AXPY_FORTRAN_IMPL(haxpydh,double,quaternion<double>,double);
  AXPY_FORTRAN_IMPL(haxpych,double,quaternion<double>,std::complex<double>);
  AXPY_FORTRAN_IMPL(haxpyhh,double,quaternion<double>,quaternion<double>);

//DOT_FORTRAN_IMPL(DOTU,hdotu,double);
  DOT_FORTRAN_IMPL(DOTC,hdotc,double);

};

