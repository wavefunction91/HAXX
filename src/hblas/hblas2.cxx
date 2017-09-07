/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#include "haxx.hpp"
#include "hblas/hblas2_impl.hpp"

namespace HAXX {

template
void HBLAS_GEMV(char, const HAXX_INT, const HAXX_INT, const double, 
  double * const, const HAXX_INT, quaternion<double> * const, const HAXX_INT, 
  const double, quaternion<double> * const, const HAXX_INT);

};

