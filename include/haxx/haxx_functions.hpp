/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HAXX_QUATERNION_FUNCTIONS_HPP
#define __INCLUDED_HAXX_QUATERNION_FUNCTIONS_HPP

#include "haxx/haxx_def.hpp"

namespace HAXX {


template <typename _F>
inline quaternion<_F> pow(const quaternion<_F>& q, const _F& a){
  // FIXME: This is wildly inefficient
  
  _F qNorm = norm(q);
  _F angle = std::acos(q.real() / qNorm);

  quaternion<_F> qIm(0.,q.imag_i(),q.imag_j(),q.imag_k());
  qIm /= (qNorm * std::sin(angle));

  quaternion<_F> qPow = std::cos(a * angle) + qIm * std::sin(a * angle);
  qPow *= std::pow(qNorm,a);

  return qPow;
};

template <typename _F>
inline quaternion<_F> sqrt(const quaternion<_F> &q) {
  return pow(q,0.5);
};

}; // namespace HAXX

#endif
