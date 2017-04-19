template <typename _F>
quaternion<_F>& quaternion<_F>::operator=(const _CF& __ct) {

  _M_real   = __ct.real();
  _M_imag_i = __ct.imag();
  _M_imag_j = _F();
  _M_imag_k = _F();

  return *this;

};

template <typename _F>
quaternion<_F>& quaternion<_F>::operator+=(const _CF& __t) {

  _M_real   += __t.real();
  _M_imag_i += __t.imag();
  return *this;

};



template <typename _F>
quaternion<_F>& quaternion<_F>::operator-=(const _CF& __t) {

  _M_real   -= __t.real();
  _M_imag_i -= __t.imag();
  return *this;

};





template <typename _F>
inline quaternion<_F> operator+(const quaternion<_F>& __x, 
  const std::complex<_F>& __y){

  quaternion<_F> __r = __x;
  __r += __y;
  return __r;

};

template <typename _F>
inline quaternion<_F> operator+(const std::complex<_F>& __x, 
  const quaternion<_F>& __y){

  quaternion<_F> __r = __y;
  __r += __x;
  return __r;

};

template <typename _F>
inline quaternion<_F> operator-(const quaternion<_F>& __x, 
  const std::complex<_F>& __y){

  quaternion<_F> __r = __x;
  __r -= __y;
  return __r;

};

template <typename _F>
inline quaternion<_F> operator-(const std::complex<_F>& __x, 
  const quaternion<_F>& __y){

  quaternion<_F> __r = -__y;
  __r += __x;
  return __r;

};



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
