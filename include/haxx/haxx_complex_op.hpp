#ifndef __INCLUDED_HAXX_COMPLEX_OP_HPP
#define __INCLUDED_HAXX_COMPLEX_OP_HPP

#include "haxx/haxx_def.hpp"

namespace HAXX {

/**
 * Assigns a complex number to this quaternion through the natural
 * embedding of the complex numbers in the quaternion algebra
 *
 * \f$ p \in \mathbb{C} \mapsto p' \in \mathbb{H} \f$ such that
 * \f$ p' = (p^R,p^I,0,0) \f$
 *
 */
template <typename _F>
quaternion<_F>& quaternion<_F>::operator=(const _CF& __ct) {

  _M_real   = __ct.real();
  _M_imag_i = __ct.imag();
  _M_imag_j = _F();
  _M_imag_k = _F();

  return *this;

};

/**
 * Adds a complex number to this quaternion in place through the
 * natural embedding of the complex numbers in the quaternion
 * algebra.
 *
 * \f$ p \in \mathbb{H}, z \in \mathbb{C} \qquad p = p + z' \f$
 * where \f$z' = (z^R, z^I,0,0) \in \mathbb{H} \f$
 *
 */
template <typename _F>
quaternion<_F>& quaternion<_F>::operator+=(const _CF& __t) {

  _M_real   += __t.real();
  _M_imag_i += __t.imag();
  return *this;

};


/**
 * Subtracts a complex number from this quaternion in place through the
 * natural embedding of the complex numbers in the quaternion
 * algebra.
 *
 * \f$ p \in \mathbb{H}, z \in \mathbb{C} \qquad p = p - z' \f$
 * where \f$z' = (z^R, z^I,0,0) \in \mathbb{H} \f$
 *
 */
template <typename _F>
quaternion<_F>& quaternion<_F>::operator-=(const _CF& __t) {

  _M_real   -= __t.real();
  _M_imag_i -= __t.imag();
  return *this;

};




/**
 *  \ingroup quaternion
 *
 *  \brief Add a complex number to a quaternion number
 *
 *  Adds a complex number and a quaternion number to return a 
 *  quaternion number through the natural embedding of the complex
 *  numbers in the quaternion algebra.
 *
 *  \f$r = p + z' \qquad 
 *    r,p\in\mathbb{H}, \quad z\in\mathbb{C}, \quad 
 *    z' = (z^R,z^I,0,0) \in\mathbb{H}\f$
 *
 */
template <typename _F>
inline quaternion<_F> operator+(const quaternion<_F>& __x, 
  const std::complex<_F>& __y){

  quaternion<_F> __r = __x;
  __r += __y;
  return __r;

};

/**
 *  \ingroup quaternion
 *
 *  \brief Add a complex number to a quaternion number
 *
 *  Adds a complex number and a quaternion number to return a 
 *  quaternion number through the natural embedding of the complex
 *  numbers in the quaternion algebra.
 *
 *  \f$r = p + z' \qquad 
 *    r,p\in\mathbb{H}, \quad z\in\mathbb{C}, \quad 
 *    z' = (z^R,z^I,0,0) \in\mathbb{H}\f$
 *
 */
template <typename _F>
inline quaternion<_F> operator+(const std::complex<_F>& __x, 
  const quaternion<_F>& __y){

  quaternion<_F> __r = __y;
  __r += __x;
  return __r;

};

/**
 *  \ingroup quaternion
 *
 *  \brief Subtract a complex number from a quaternion number
 *
 *  Subtracts a complex number from a quaternion number to return a 
 *  quaternion number through the natural embedding of the complex
 *  numbers in the quaternion algebra.
 *
 *  \f$r = p - z' \qquad 
 *    r,p\in\mathbb{H}, \quad z\in\mathbb{C}, \quad 
 *    z' = (z^R,z^I,0,0) \in\mathbb{H}\f$
 *
 */
template <typename _F>
inline quaternion<_F> operator-(const quaternion<_F>& __x, 
  const std::complex<_F>& __y){

  quaternion<_F> __r = __x;
  __r -= __y;
  return __r;

};

/**
 *  \ingroup quaternion
 *
 *  \brief Subtract a quaternion number from a complex number
 *
 *  Subtracts a quaternion number from a complex number to return a 
 *  quaternion number through the natural embedding of the complex
 *  numbers in the quaternion algebra.
 *
 *  \f$r = z' - p \qquad 
 *    r,p\in\mathbb{H}, \quad z\in\mathbb{C}, \quad 
 *    z' = (z^R,z^I,0,0) \in\mathbb{H}\f$
 *
 */
template <typename _F>
inline quaternion<_F> operator-(const std::complex<_F>& __x, 
  const quaternion<_F>& __y){

  quaternion<_F> __r = -__y;
  __r += __x;
  return __r;

};



/**
 *  \ingroup quaternion
 *
 *  \brief Right multiply a quaternion number by a complex number
 *
 *  Right multiplies a quaternion number by a complex number to return a 
 *  quaternion number through the natural embedding of the complex
 *  numbers in the quaternion algebra. 
 *
 *  \f$r = pz' \qquad 
 *    r,p\in\mathbb{H}, \quad z\in\mathbb{C}, \quad 
 *    z' = (z^R,z^I,0,0) \in\mathbb{H}\f$
 *
 */
template <typename _F>
inline quaternion<_F> operator*(const quaternion<_F>& __x,
  const std::complex<_F>& __y) {

  quaternion<_F> __r;

  // This is a really naive algorithm
  __r.real(__x.real()   * __y.real()   - __x.imag_i() * __y.imag());
  __r.imag_i(__x.real()   * __y.imag() + __x.imag_i() * __y.real());
  __r.imag_j(__x.imag_j() * __y.real() + __x.imag_k() * __y.imag());
  __r.imag_k(__x.imag_k() * __y.real() - __x.imag_j() * __y.imag());
  
  return __r;
};

/**
 *  \ingroup quaternion
 *
 *  \brief Left multiply a quaternion number by a complex number
 *
 *  Left multiplies a quaternion number by a complex number to return a 
 *  quaternion number through the natural embedding of the complex
 *  numbers in the quaternion algebra. 
 *
 *  \f$r = pz' \qquad 
 *    r,p\in\mathbb{H}, \quad z\in\mathbb{C}, \quad 
 *    z' = (z^R,z^I,0,0) \in\mathbb{H}\f$
 *
 */
template <typename _F>
inline quaternion<_F> operator*(const std::complex<_F>& __x,
  const quaternion<_F>& __y) {

  quaternion<_F> __r;

  // This is a really naive algorithm
  __r.real(__x.real()   * __y.real()   - __x.imag() * __y.imag_i());
  __r.imag_i(__x.real() * __y.imag_i() + __x.imag() * __y.real());
  __r.imag_j(__x.real() * __y.imag_j() - __x.imag() * __y.imag_k());
  __r.imag_k(__x.real() * __y.imag_k() + __x.imag() * __y.imag_j());
  
  return __r;
};


}; // HAXX namespace

#endif
