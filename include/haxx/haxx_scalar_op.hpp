#ifndef __INCLUDED_HAXX_SCALAR_OP_HPP
#define __INCLUDED_HAXX_SCALAR_OP_HPP

#include "haxx/haxx_def.hpp"

namespace HAXX {

/**
 * Assigns a real number to this quaternion through the natural
 * embedding of the real numbers in the quaternion algebra
 *
 * \f$ p \in \mathbb{R} \mapsto p' \in \mathbb{H} \f$ such that
 * \f$ p' = (p,0,0,0) \f$
 *
 */
template <typename _F>
quaternion<_F>& quaternion<_F>::operator=(const _F& __t) {

  _M_real   = __t;
  _M_imag_i = _F();
  _M_imag_j = _F();
  _M_imag_k = _F();

  return *this;

};



/**
 * Adds a real number to this quaternion in place through the
 * natural embedding of the real numbers in the quaternion
 * algebra.
 *
 * \f$ p \in \mathbb{H}, a \in \mathbb{R} \qquad p = p + a' \f$
 * where \f$a' = (a, 0,0,0) \in \mathbb{H} \f$
 *
 */
template <typename _F>
quaternion<_F>& quaternion<_F>::operator+=(const _F& __t) {

  _M_real += __t;
  return *this;

};

/**
 * Subtracts a real number from this quaternion in place through the
 * natural embedding of the real numbers in the quaternion
 * algebra.
 *
 * \f$ p \in \mathbb{H}, a \in \mathbb{R} \qquad p = p - a' \f$
 * where \f$a' = (a,0,0,0) \in \mathbb{H} \f$
 *
 */
template <typename _F>
quaternion<_F>& quaternion<_F>::operator-=(const _F& __t) {

  _M_real -= __t;
  return *this;

};

/**
 * Multiplies this quaternion by a scalar in place through the
 * natural embedding of the real numbers in the quaternion
 * algebra.
 *
 * \f$ p \in \mathbb{H}, a \in \mathbb{R} \qquad p = pa' =(ap^R,ap^I,ap^J,ap^K)\f$
 * where \f$a' = (a,0,0,0) \in \mathbb{H} \f$
 *
 */
template <typename _F>
quaternion<_F>& quaternion<_F>::operator*=(const _F& __t) {

  _M_real   *= __t;
  _M_imag_i *= __t;
  _M_imag_j *= __t;
  _M_imag_k *= __t;

  return *this;

};

/**
 * Divides this quaternion by a scalar in place through the
 * natural embedding of the real numbers in the quaternion
 * algebra.
 *
 * \f$ p \in \mathbb{H}, a \in \mathbb{R} \qquad 
 *  p = p(a')^{-1} =(p^R/a,p^I/a,p^J/a,p^K/a)\f$
 * where \f$a' = (a,0,0,0) \in \mathbb{H} \f$
 *
 */
template <typename _F>
quaternion<_F>& quaternion<_F>::operator/=(const _F& __t) {

  _M_real   /= __t;
  _M_imag_i /= __t;
  _M_imag_j /= __t;
  _M_imag_k /= __t;

  return *this;

};





/**
 *  \ingroup quaternion
 *
 *  \brief Add a real number to a quaternion number
 *
 *  Adds a real number and a quaternion number to return a 
 *  quaternion number through the natural embedding of the real
 *  numbers in the quaternion algebra.
 *
 *  \f$r = p + a' \qquad 
 *    r,p\in\mathbb{H}, \quad a\in\mathbb{R}, \quad 
 *    a' = (a,0,0,0) \in\mathbb{H}\f$
 *
 */
template <typename _F>
inline quaternion<_F> operator+(const quaternion<_F>& __x, const _F& __y){

  quaternion<_F> __r = __x;
  __r += __y;
  return __r;

};

/**
 *  \ingroup quaternion
 *
 *  \brief Add a real number to a quaternion number
 *
 *  Adds a real number and a quaternion number to return a 
 *  quaternion number through the natural embedding of the real
 *  numbers in the quaternion algebra.
 *
 *  \f$r = a' + p \qquad 
 *    r,p\in\mathbb{H}, \quad a\in\mathbb{R}, \quad 
 *    a' = (a,0,0,0) \in\mathbb{H}\f$
 *
 */
template <typename _F>
inline quaternion<_F> operator+(const _F& __x, const quaternion<_F>& __y){

  quaternion<_F> __r = __y;
  __r += __x;
  return __r;

};

/**
 *  \ingroup quaternion
 *
 *  \brief Subtract a real number from a quaternion number
 *
 *  Subtract a real number from a quaternion number to return a 
 *  quaternion number through the natural embedding of the real
 *  numbers in the quaternion algebra.
 *
 *  \f$r = p - a' \qquad 
 *    r,p\in\mathbb{H}, \quad a\in\mathbb{R}, \quad 
 *    a' = (a,0,0,0) \in\mathbb{H}\f$
 *
 */
template <typename _F>
inline quaternion<_F> operator-(const quaternion<_F>& __x, const _F& __y){

  quaternion<_F> __r = __x;
  __r -= __y;
  return __r;

};

/**
 *  \ingroup quaternion
 *
 *  \brief Subtract a quaternion number from a real number
 *
 *  Subtract a quaternion number from a real number to return a 
 *  quaternion number through the natural embedding of the real
 *  numbers in the quaternion algebra.
 *
 *  \f$r = a' - p \qquad 
 *    r,p\in\mathbb{H}, \quad a\in\mathbb{R}, \quad 
 *    a' = (a,0,0,0) \in\mathbb{H}\f$
 *
 */
template <typename _F>
inline quaternion<_F> operator-(const _F& __x, const quaternion<_F>& __y){

  quaternion<_F> __r = -__y;
  __r += __x;
  return __r;

};


/**
 * \ingroup quaternion
 *
 * \brief Multply a quaternion and a scalar
 *
 * Multiplies a quaternion by a scalar through the
 * natural embedding of the real numbers in the quaternion
 * algebra.
 *
 * \f$ r,p \in \mathbb{H}, a \in \mathbb{R} \qquad r = pa' =(ap^R,ap^I,ap^J,ap^K)\f$
 * where \f$a' = (a,0,0,0) \in \mathbb{H} \f$
 *
 */
template <typename _F>
inline quaternion<_F> operator*(const quaternion<_F>& __x, const _F& __y){

  quaternion<_F> __r = __x;
  __r *= __y;
  return __r;

};

/**
 * \ingroup quaternion
 *
 * \brief Multply a quaternion and a scalar
 *
 * Multiplies a quaternion by a scalar through the
 * natural embedding of the real numbers in the quaternion
 * algebra.
 *
 * \f$ p,r \in \mathbb{H}, a \in \mathbb{R} \qquad r = a'p =(ap^R,ap^I,ap^J,ap^K)\f$
 * where \f$a' = (a,0,0,0) \in \mathbb{H} \f$
 *
 */
template <typename _F>
inline quaternion<_F> operator*(const _F& __x, const quaternion<_F>& __y){

  quaternion<_F> __r = __y;
  __r *= __x;
  return __r;

};

/**
 * \ingroup quaternion
 *
 * \brief Divide a quaternion by a scalar
 *
 * Divides a quaternion by a scalar through the
 * natural embedding of the real numbers in the quaternion
 * algebra.
 *
 * \f$ p,r \in \mathbb{H}, a \in \mathbb{R} \qquad r = p(a')^{-1} =(p^R/a,p^I/a,p^J/a,p^K/a)\f$
 * where \f$a' = (a,0,0,0) \in \mathbb{H} \f$
 *
 */
template <typename _F>
inline quaternion<_F> operator/(const quaternion<_F>& __x, const _F& __y){

  quaternion<_F> __r = __x;
  __r /= __y;
  return __r;

};

/**
 * \ingroup quaternion
 *
 * \brief Divide a scalar by a quaternion
 *
 * Divides a scalar by a quaternion through the
 * natural embedding of the real numbers in the quaternion
 * algebra.
 *
 * \f$ p,r \in \mathbb{H}, a \in \mathbb{R} \qquad r = a'p^{-1} \f$
 * where \f$a' = (a,0,0,0) \in \mathbb{H} \f$
 *
 */
template <typename _F>
inline quaternion<_F> operator/(const _F& __x, const quaternion<_F>& __y){

  quaternion<_F> __r = inv(__y);
  __r *= __x;
  return __r;

};

}; // HAXX namespace

#endif
