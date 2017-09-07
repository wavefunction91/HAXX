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

namespace HAXX {

template
void HBLAS_SWAP(HAXX_INT, quaternion<double> *, HAXX_INT, 
  quaternion<double> *, HAXX_INT);

template
void HBLAS_COPY(HAXX_INT, quaternion<double> *, HAXX_INT, 
  quaternion<double> *, HAXX_INT);

};

