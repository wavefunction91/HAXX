#ifndef __INCLUDED_HAXX_SCALAR_OP_HPP
#define __INCLUDED_HAXX_SCALAR_OP_HPP

#include "haxx/haxx_def.hpp"

namespace HAXX {

template <typename _F>
quaternion<_F>& quaternion<_F>::operator=(const _F& __t) {

  _M_real   = __t;
  _M_imag_i = _F();
  _M_imag_j = _F();
  _M_imag_k = _F();

  return *this;

};



template <typename _F>
quaternion<_F>& quaternion<_F>::operator+=(const _F& __t) {

  _M_real += __t;
  return *this;

};

template <typename _F>
quaternion<_F>& quaternion<_F>::operator-=(const _F& __t) {

  _M_real -= __t;
  return *this;

};

template <typename _F>
quaternion<_F>& quaternion<_F>::operator*=(const _F& __t) {

  _M_real   *= __t;
  _M_imag_i *= __t;
  _M_imag_j *= __t;
  _M_imag_k *= __t;

  return *this;

};

template <typename _F>
quaternion<_F>& quaternion<_F>::operator/=(const _F& __t) {

  _M_real   /= __t;
  _M_imag_i /= __t;
  _M_imag_j /= __t;
  _M_imag_k /= __t;

  return *this;

};





template <typename _F>
inline quaternion<_F> operator+(const quaternion<_F>& __x, const _F& __y){

  quaternion<_F> __r = __x;
  __r += __y;
  return __r;

};

template <typename _F>
inline quaternion<_F> operator+(const _F& __x, const quaternion<_F>& __y){

  quaternion<_F> __r = __y;
  __r += __x;
  return __r;

};

template <typename _F>
inline quaternion<_F> operator-(const quaternion<_F>& __x, const _F& __y){

  quaternion<_F> __r = __x;
  __r -= __y;
  return __r;

};

template <typename _F>
inline quaternion<_F> operator-(const _F& __x, const quaternion<_F>& __y){

  quaternion<_F> __r = -__y;
  __r += __x;
  return __r;

};

template <typename _F>
inline quaternion<_F> operator*(const quaternion<_F>& __x, const _F& __y){

  quaternion<_F> __r = __x;
  __r *= __y;
  return __r;

};

template <typename _F>
inline quaternion<_F> operator*(const _F& __x, const quaternion<_F>& __y){

  quaternion<_F> __r = __y;
  __r *= __x;
  return __r;

};

template <typename _F>
inline quaternion<_F> operator/(const quaternion<_F>& __x, const _F& __y){

  quaternion<_F> __r = __x;
  __r /= __y;
  return __r;

};

template <typename _F>
inline quaternion<_F> operator/(const _F& __x, const quaternion<_F>& __y){

  quaternion<_F> __r = inv(__y);
  __r *= __x;
  return __r;

};

}; // HAXX namespace

#endif
