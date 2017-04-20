#ifndef __INCLUDED_HAXX_QUATERNION_REDUCTION_HPP
#define __INCLUDED_HAXX_QUATERNION_REDUCTION_HPP

#include "haxx/haxx_def.hpp"

namespace HAXX {

/**
 *  \ingroup quaternion
 *
 *  \brief Returns the norm of a quaternion
 *
 *  \f$ \vert\vert q \vert\vert = \sqrt{a^2 + b^2 + c^2 + d^2} \qquad q = (a,b,c,d) \in \mathbb{H}\f$
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
 *  \ingroup quaternion
 *
 *  \brief Returns quaternion conjugate of a quaternion
 *
 *  \f$ q^* = (q^R,-q^I,-q^J,-q^K)\f$
 *
 */
template <typename _F>
inline quaternion<_F> conj(const quaternion<_F>& __q) {

  return quaternion<_F>(__q.real(),-__q.imag_i(),-__q.imag_j(),-__q.imag_k());

};

/**
 *  \ingroup quaternion
 *
 *  \brief Return the inverse of a quaternion
 *
 *  \f$ q^{-1} = \dfrac{q^*}{\vert\vert q \vert\vert^2} \f$
 *
 */
template <typename _F>
inline quaternion<_F> inv(const quaternion<_F>& __q) {

  _F nrm = norm(__q);
  return conj(__q) / nrm / nrm;

};

};

#endif
