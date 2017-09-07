/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HAXX_QUATERNION_REDUCTION_HPP
#define __INCLUDED_HAXX_QUATERNION_REDUCTION_HPP

#include "haxx/haxx_def.hpp"

namespace HAXX {

/**
 *  \f$ \vert\vert q \vert\vert = \sqrt{a^2 + b^2 + c^2 + d^2} 
 *    \qquad q = (a,b,c,d) \in \mathbb{H}\f$
 *
 *  Return type is the value_type
 */
template <typename _F>
inline _F norm(const quaternion<_F>& __q) {
  _F nmsq = __q.real() * __q.real();
  nmsq += __q.imag_i() * __q.imag_i();
  nmsq += __q.imag_j() * __q.imag_j();
  nmsq += __q.imag_k() * __q.imag_k();

  return std::sqrt(nmsq);
};

/**
 *  \f$ q^* = (q^R,-q^I,-q^J,-q^K)\f$
 */
template <typename _F>
inline quaternion<_F> conj(const quaternion<_F>& __q) {

  return quaternion<_F>(__q.real(),-__q.imag_i(),-__q.imag_j(),-__q.imag_k());

};

/**
 *  \f$ q^{-1} = \dfrac{q^*}{\vert\vert q \vert\vert^2} \f$
 */
template <typename _F>
inline quaternion<_F> inv(const quaternion<_F>& __q) {

  _F nrm = norm(__q);
  return conj(__q) / nrm / nrm;

};

/**
 *  \f$ r = \dfrac{q}{\vert\vert q \vert\vert} \f$
 */
template <typename _F>
inline quaternion<_F> versor(const quaternion<_F>& __q) {

  _F nrm = norm(__q);
  return __q / nrm;

};

/**
 *  \f$ [p,q] = pq - qp \f$
 */
template <typename _F>
inline quaternion<_F> comm(const quaternion<_F>& p, const quaternion<_F>& q) {

  return p * q - q * p;

};

template <typename _F> 
inline quaternion<_F> comm(const quaternion<_F>& p, const std::complex<_F>& q) {
/*
  std::complex<double> jPt(p.imag_j(),p.imag_k());
  return quaternion<_F>(std::complex<double>(0.),2.*jPt*std::imag(q));
*/
  return p * q - q * p;
};

template <typename _F> 
inline quaternion<_F> comm(const std::complex<_F>& p, const quaternion<_F>& q) {

  return -comm(q,p);

}


template <typename _F> 
inline quaternion<_F> comm(const _F& p, const quaternion<_F>& q) {

  return quaternion<_F>(0.);

}

template <typename _F> 
inline quaternion<_F> comm(const quaternion<_F>& q, const _F& p) {

  return quaternion<_F>(0.);

}


};

#endif
