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

};

#endif
