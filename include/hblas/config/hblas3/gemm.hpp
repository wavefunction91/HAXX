/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HBLAS_CONFIG_HBLAS3_GEMM_HPP
#define __INCLUDED_HBLAS_CONFIG_HBLAS3_GEMM_HPP


#include "hblas/config/types.hpp"

// Caching dimensions
#define MC 64
#define NC 1024
#define KC 64

// Register block size
#define MR 2
#define NR 2


// Determine where to factor ALPHA scaling
#define _FACTOR_ALPHA_IN_A_PACK
//#define _FACTOR_ALPHA_IN_B_PACK

#if defined(_FACTOR_ALPHA_IN_A_PACK) && defined(_FACTOR_ALPHA_IN_B_PACK)
  #error "Cannot factor ALPHA into both A and B packs"
#endif


// Determine where to factor the transpose operation
// for GEMM kernel
#define _FACTOR_TRANSPOSE_INTO_A_PACK
#define _FACTOR_TRANSPOSE_INTO_B_PACK



// Determine packing utility for GEMM

#include "hblas/pack/pack.hpp"

#ifdef _FACTOR_TRANSPOSE_INTO_B_PACK

  #define BPACKT  NPACK< NR, _AMATF, GenericPackOps_T2<> >    
  #define BPACKCT NPACK< NR, _AMATF, ConjPackOps_T2   <> >
  #define BPACKR  TPACK< NR, _AMATF, ConjPackOps_T2   <> >
  #define BPACK   TPACK< NR, _AMATF, GenericPackOps_T2<> >      

#else

  #define BPACKT  NPACK< NR, _BMATF, GenericPackOps<_BMATF> >
  #define BPACKCT NPACK< NR, _BMATF, ConjPackOps   <_BMATF> >
  #define BPACKR  TPACK< NR, _BMATF, ConjPackOps   <_BMATF> >
  #define BPACK   TPACK< NR, _BMATF, GenericPackOps<_BMATF> >   

#endif

#ifdef _FACTOR_TRANSPOSE_INTO_A_PACK

  #define APACKT  TPACK< MR, _AMATF, GenericPackOps_T1<> > 
  #define APACKCT TPACK< MR, _AMATF, ConjPackOps_T1   <> >
  #define APACKR  NPACK< MR, _AMATF, ConjPackOps_T1   <> >
  #define APACK   NPACK< MR, _AMATF, GenericPackOps_T1<> >   

#else

  #define APACKT  TPACK< MR, _AMATF, GenericPackOps<_AMATF> >
  #define APACKCT TPACK< MR, _AMATF, ConjPackOps   <_AMATF> >
  #define APACKR  NPACK< MR, _AMATF, ConjPackOps   <_AMATF> >
  #define APACK   NPACK< MR, _AMATF, GenericPackOps<_AMATF> >    

#endif

#endif
