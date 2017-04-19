#ifndef __INCLUDED_HAXX_DEF_HPP
#define __INCLUDED_HAXX_DEF_HPP


#include <complex> // Support for complex numbers

// HAXX Namespace
namespace HAXX {

// Primary template class quaternion
template <typename _F>
struct quaternion {

  /// Value typedef
  typedef _F value_type;
  typedef std::complex<_F> _CF;

  /// Default Constructor.
  /// Constructs a quaternion from real quadruple
  ///   q = __r + __i * i + __j * j + __k * k
  constexpr quaternion(const _F& __r = _F(), const _F& __i = _F(), 
    const _F& __j = _F(), const _F& __k = _F()) :
    _M_real(__r), _M_imag_i(__i), _M_imag_j(__j), _M_imag_k(__k){ } 

  /// Complex Constructor.
  /// Constructs a quaternion using Cayley-Dickson construction of
  ///   the complex numbers
  ///
  ///   q = __a + __b * j
  constexpr quaternion(const _CF& __a, const _CF& __b = _CF()) :
    _M_real(__a.real()), _M_imag_i(__a.imag()),
    _M_imag_j(__b.real()), _M_imag_k(__b.imag()){ }

  /// Copy Constructor.
  constexpr quaternion(const quaternion&) = default;


  /// Converting Constructor.
  template <typename _G>
    constexpr quaternion(const quaternion<_G>& __q) :
      _M_real(__q.real()), _M_imag_i(__q.imag_i()),
      _M_imag_j(__q.imag_j()), _M_imag_k(__q.imag_k()){ }  






  // Getters for components of the quaternion
  constexpr _F real()   const { return _M_real;   }
  constexpr _F imag_i() const { return _M_imag_i; }
  constexpr _F imag_j() const { return _M_imag_j; }
  constexpr _F imag_k() const { return _M_imag_k; }


  // Setters for components of the quaternion
  void real(_F __x)   { _M_real   = __x; }
  void imag_i(_F __x) { _M_imag_i = __x; }
  void imag_j(_F __x) { _M_imag_j = __x; }
  void imag_k(_F __x) { _M_imag_k = __x; }

  // Real / Complex Assignment Operators

  /// Assign a scalar to this quaternion number
  quaternion<_F>& operator=(const _F&);

  /// Assign a complex number to this quaternion number
  quaternion<_F>& operator=(const _CF&);





  // Real / Complex Arithmetatic operators

  /// Add a scalar to this quaternion number
  quaternion<_F>& operator+=(const _F&);

  /// Subtract a scalar from this quaternion number
  quaternion<_F>& operator-=(const _F&);

  /// Multiply a scalar by this quaternion number
  quaternion<_F>& operator*=(const _F&);

  /// Divide this quaternion number by a scalar
  quaternion<_F>& operator/=(const _F&);


  /// Add a complex number to this quaternion number
  quaternion<_F>& operator+=(const _CF&);

  /// Subtract a complex number from this quaternion number
  quaternion<_F>& operator-=(const _CF&);

  // Self multiplication by a complex number is ambiguous







  // Quaternion operators

  /// Assign a quaternion number (same type) to this quaternion number
  quaternion& operator=(const quaternion&) = default;

  /// Assign a quaternion number (different type) to this quaternion number
  template <typename _G> quaternion<_F>& operator=(const quaternion<_G>&);

  /// Add a quaternion number to this quaternion number
  template <typename _G> quaternion<_F>& operator+=(const quaternion<_G>&);

  /// Subtract a quaternion number from this quaternion number
  template <typename _G> quaternion<_F>& operator-=(const quaternion<_G>&);

  // Self multiplication by a quaternion number is ambiguous






  constexpr quaternion __rep() const { return *this; }

  private:
    _F _M_real;   /// Real (Scalar) part of the quaternion
    _F _M_imag_i; /// Imaginary (i) part of the quaternion
    _F _M_imag_j; /// Imaginary (j) part of the quaternion
    _F _M_imag_k; /// Imaginary (k) part of the quaternion



}; // struct quaternion


}; // HAXX namespace


#endif
