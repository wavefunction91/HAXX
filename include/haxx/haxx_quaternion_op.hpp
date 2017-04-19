
namespace HAXX {

template <typename _F>
template <typename _G>
quaternion<_F>& quaternion<_F>::operator=(const quaternion<_G> &__q) {

  _M_real   = __q.real();
  _M_imag_i = __q.imag_i();
  _M_imag_j = __q.imag_j();
  _M_imag_k = __q.imag_k();

  return *this;

};

template <typename _F>
template <typename _G>
quaternion<_F>& quaternion<_F>::operator+=(const quaternion<_G> &__q) {

  _M_real   += __q.real();
  _M_imag_i += __q.imag_i();
  _M_imag_j += __q.imag_j();
  _M_imag_k += __q.imag_k();

  return *this;

};

template <typename _F>
template <typename _G>
quaternion<_F>& quaternion<_F>::operator-=(const quaternion<_G> &__q) {

  _M_real   -= __q.real();
  _M_imag_i -= __q.imag_i();
  _M_imag_j -= __q.imag_j();
  _M_imag_k -= __q.imag_k();

  return *this;

};




template <typename _F>
inline quaternion<_F> operator-(const quaternion<_F>& __x) {

  return quaternion<_F>(-__x.real(),-__x.imag_i(),-__x.imag_j(),-__x.imag_k());

}






template <typename _F>
inline quaternion<_F> operator+(const quaternion<_F>& __x, 
  const quaternion<_F>& __y) {

  quaternion<_F> __r = __x;
  __r += __y;
  return __r;

};

template <typename _F>
inline quaternion<_F> operator-(const quaternion<_F>& __x, 
  const quaternion<_F>& __y) {

  quaternion<_F> __r = __x;
  __r -= __y;
  return __r;

};

template <typename _F>
inline quaternion<_F> operator*(const quaternion<_F>& __x,
  const quaternion<_F>& __y) {

  quaternion<_F> __r;

  // This is a really naive algorithm
  __r.real(__x.real()   * __y.real()   - __x.imag_i() * __y.imag_i() -
           __x.imag_j() * __y.imag_j() - __x.imag_k() * __y.imag_k());

  __r.imag_i(__x.real()   * __y.imag_i() + __x.imag_i() * __y.real()   +
             __x.imag_j() * __y.imag_k() - __x.imag_k() * __y.imag_j());

  __r.imag_j(__x.real()   * __y.imag_j() - __x.imag_i() * __y.imag_k()  +
             __x.imag_j() * __y.real()   + __x.imag_k() * __y.imag_i());

  __r.imag_k(__x.real()   * __y.imag_k() + __x.imag_i() * __y.imag_j()  -
             __x.imag_j() * __y.imag_i() + __x.imag_k() * __y.real());
  
  return __r;
};


}; // HAXX namespace
