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
void HBLAS_GEMV(char, HAXX_INT, HAXX_INT, double, double *, 
  HAXX_INT, quaternion<double> *, HAXX_INT, double, 
  quaternion<double> *, HAXX_INT);

};

