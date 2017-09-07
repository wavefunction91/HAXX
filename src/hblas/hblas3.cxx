/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#include "haxx.hpp"
#include "hblas/hblas3_impl.hpp"

namespace HAXX {


// Instantiate REAL-QUATERION matrix multiplication from template
template
void HBLAS_GEMM(const char, const char, const HAXX_INT, const HAXX_INT, 
  const HAXX_INT, const double, double * const, const HAXX_INT, 
  quaternion<double> * const, const HAXX_INT, const double, 
  quaternion<double> * const, const HAXX_INT);

};

