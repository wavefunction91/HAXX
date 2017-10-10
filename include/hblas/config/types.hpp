/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_COMPILE_CONFIG_HPP
#define __INCLUDED_HBLAS_COMPILE_CONFIG_HPP

#include "util/constants.hpp"

// Determine type of scaling parameter ALPHA
#ifdef ALPHAF
  #if ALPHAF == DOUBLE
    #define _ALPHAF double
  #elif ALPHAF == DCOMPLEX
    #define _ALPHAF std::complex<double>
  #elif ALPHAF == DQUATERNION
    #define _ALPHAF quaternion<double>
  #endif
#endif

// Determine type of scaling parameter BETA
#ifdef BETAF
  #if BETAF == DOUBLE
    #define _BETAF double
  #elif BETAF == DCOMPLEX
    #define _BETAF std::complex<double>
  #elif BETAF == DQUATERNION
    #define _BETAF quaternion<double>
  #endif
#endif

// Determine type of matrix type AMAT
#ifdef AMATF
  #if AMATF == DOUBLE
    #define _AMATF double
  #elif AMATF == DCOMPLEX
    #define _AMATF std::complex<double>
  #elif AMATF == DQUATERNION
    #define _AMATF quaternion<double>
  #endif
#endif

// Determine type of matrix type BMAT
#ifdef BMATF
  #if BMATF == DOUBLE
    #define _BMATF double
  #elif BMATF == DCOMPLEX
    #define _BMATF std::complex<double>
  #elif BMATF == DQUATERNION
    #define _BMATF quaternion<double>
  #endif
#endif



#endif
