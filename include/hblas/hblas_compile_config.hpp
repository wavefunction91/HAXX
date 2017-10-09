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

/** COMMON DEFS FOR ALL HBLAS ROUTINES **/

// Determine type of scaling parameter ALPHA
#ifdef ALPHAF
  #if ALPHAF == DOUBLE
    #define _ALPHAF double
  #elif ALPHAF == DCOMPLEX
    #define _ALPHAF std::complex<double>
  #elif ALPHAF == DQUATERNION
    #define _ALPHAF quaternion<double>
  #else
    #error GEMM Only Supports 64-bit floats
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
  #else
    #error GEMM Only Supports 64-bit floats
  #endif
#endif



#endif
