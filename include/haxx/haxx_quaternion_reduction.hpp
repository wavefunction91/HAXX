
namespace HAXX {

template <typename _F>
inline _F norm(const quaternion<_F>& __q) {
  _F nmsq = __q.real() * __q.real();
  nmsq += __q.imag_i() * __q.imag_i();
  nmsq += __q.imag_j() * __q.imag_j();
  nmsq += __q.imag_k() * __q.imag_k();

  return std::sqrt(nmsq);
};

template <typename _F>
inline quaternion<_F> conj(const quaternion<_F>& __q) {

  return quaternion<_F>(__q.real(),-__q.imag_i(),-__q.imag_j(),-__q.imag_k());

};


template <typename _F>
inline quaternion<_F> inv(const quaternion<_F>& __q) {

  _F nrm = norm(__q);
  return conj(__q) / nrm / nrm;

};

};
