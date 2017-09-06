/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_HAXX_DEF_HPP
#define __INCLUDED_HAXX_DEF_HPP


#include <complex> // Support for complex numbers
#include <string>
#include <sstream>
#include <array>

// HAXX Namespace
namespace HAXX {

/**
 *  @defgroup quaternion The Quaternion Struct
 *  The templated quaternion struct definition and associated
 *  arithmetic functions for the quaternion algebra.
 *
 *  @{
 */

/**
 * Primary template class quaternion
 *
 */
template <typename _F>
struct quaternion {

  /// Value typedef
  typedef _F value_type;

  /// Complex value typedef
  typedef std::complex<_F> _CF;

  /// \brief Default Constructor.
  /// Constructs a quaternion from real quadruple
  ///
  ///   \f$q = a + bi + cj + dk = (a,b,c,d) \qquad a,b,c,d\in\mathbb{R}\f$
  constexpr quaternion(const _F& a = _F(), const _F& b = _F(), 
    const _F& c = _F(), const _F& d = _F()) :
    _M_real(a), _M_imag_i(b), _M_imag_j(c), _M_imag_k(d){ } 

  /// \brief Complex Constructor.
  /// Constructs a quaternion using Cayley-Dickson construction of
  ///   the complex numbers
  ///
  ///   \f$ q = a + bj = a^R + a^Ii + b^Rj + b^Ik \qquad a,b \in \mathbb{C}\f$
  constexpr quaternion(const _CF& a, const _CF& b = _CF()) :
    _M_real(a.real()), _M_imag_i(a.imag()),
    _M_imag_j(b.real()), _M_imag_k(b.imag()){ }


  /// \brief Real Array Constructor
  /// Constructs a quaternion from a std::array<_F,4>
  constexpr quaternion(const std::array<_F,4>& q) :
    quaternion(q[0],q[1],q[2],q[3]){ };

  /// \brief Complex Array Constructor
  /// Constructs a quaternion from a std::array<_CF,4>
  constexpr quaternion(const std::array<_CF,2>& q) :
    quaternion(q[0],q[1]){ };

  /// Copy Constructor.
  constexpr quaternion(const quaternion&) = default;

  ~quaternion() = default;

  /// Converting Constructor.
  template <typename _G>
    constexpr quaternion(const quaternion<_G>& __q) :
      _M_real(__q.real()), _M_imag_i(__q.imag_i()),
      _M_imag_j(__q.imag_j()), _M_imag_k(__q.imag_k()){ }  






  // Getters for components of the quaternion
 
  /// Returns the real (scalar) part of the quaternion
  constexpr _F real()   const { return _M_real;   }
  /// Returns the imaginary (i) part of the quaternion
  constexpr _F imag_i() const { return _M_imag_i; }
  /// Returns the imaginary (j) part of the quaternion
  constexpr _F imag_j() const { return _M_imag_j; }
  /// Returns the imaginary (k) part of the quaternion
  constexpr _F imag_k() const { return _M_imag_k; }
  /// Returns the imaginary vector {i,j,k} of the quaternion
  constexpr std::array<_F,3> imag() const { return { _M_imag_i,_M_imag_j,_M_imag_k }; }


  // Setters for components of the quaternion

  /// Sets the real (scalar) part of the quaternion
  void real(_F __x)   { _M_real   = __x; }
  /// Sets the imaginary (i) part of the quaternion
  void imag_i(_F __x) { _M_imag_i = __x; }
  /// Sets the imaginary (j) part of the quaternion
  void imag_j(_F __x) { _M_imag_j = __x; }
  /// Sets the imaginary (k) part of the quaternion
  void imag_k(_F __x) { _M_imag_k = __x; }
  /// Sets the imaginary {i,j,k} part of the quaternion
  void imag(const std::array<_F,3>& __x) {
    _M_imag_i = __x[0];
    _M_imag_j = __x[1];
    _M_imag_k = __x[2];
  }

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

  /// Quaternion conversion assignment
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

///  Insertion operator for quaternion values.
template<typename _F, typename _CharT, class _Traits>
  std::basic_ostream<_CharT, _Traits>&
  operator<<(std::basic_ostream<_CharT, _Traits>& __os, const quaternion<_F>& __x)
  {
    std::basic_ostringstream<_CharT, _Traits> __s;
    __s.flags(__os.flags());
    __s.imbue(__os.getloc());
    __s.precision(__os.precision());
    __s << '(' << __x.real() << ',' << __x.imag_i() 
        << ',' << __x.imag_j() << ',' << __x.imag_k() << ')';
    return __os << __s.str();
  }

template <typename _F>
std::string to_string(quaternion<_F> &x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
};


/**
 *  @defgoup QOp Quaterion Unary Operations
 *  Unary operations on quaternion numbers
 *  @{
 */

/// Extract real part of a quaternion number
template <typename _F>
  inline _F real(const quaternion<_F> &x){ return x.real(); } 

/// Extract i-imaginary part of a quaternion number
template <typename _F>
  inline _F imag_i(const quaternion<_F> &x){ return x.imag_i(); } 

/// Extract j-imaginary part of a quaternion number
template <typename _F>
  inline _F imag_j(const quaternion<_F> &x){ return x.imag_j(); } 

/// Extract k-imaginary part of a quaternion number
template <typename _F>
  inline _F imag_k(const quaternion<_F> &x){ return x.imag_k(); } 


/* @} */ // QOp

/**
 *  @defgroup QQOp Quaternion - Quaternion Binary Operations
 *  Binary arithmetic operations between quaternion numbers
 *  @{
 */

/// Add two quaternion numbers
template <typename _F>
inline quaternion<_F> operator+(const quaternion<_F>&, const quaternion<_F>&);

/// Subtract two quaternion numbers
template <typename _F>
inline quaternion<_F> operator-(const quaternion<_F>&, const quaternion<_F>&);

/// Multiply two quaternion numbers
template <typename _F>
inline quaternion<_F> operator*(const quaternion<_F>&, const quaternion<_F>&);

/// Negate a quaternion number
template <typename _F>
inline quaternion<_F> operator-(const quaternion<_F>& __x);


/// Check the equality of two quaternions
template <typename _F>
inline bool operator==(const quaternion<_F>&, const quaternion<_F>&);
/// Check the inequality of two quaternions
template <typename _F>
inline bool operator==(const quaternion<_F>&, const quaternion<_F>&);

/* @} */ // QQOp


/**
 *  @defgroup QROp Quaternion - Real Binary Operations
 *  Binary arithmetic operations between quaternion and real numbers
 *  @{
 */

/// Add a real number to a quaternion number
template <typename _F>
inline quaternion<_F> operator+(const quaternion<_F>&, const _F&);

/// Add a real number to a quaternion number
template <typename _F>
inline quaternion<_F> operator+(const _F&, const quaternion<_F>&);

/// Subtract a real number from a quaternion number
template <typename _F>
inline quaternion<_F> operator-(const quaternion<_F>&, const _F&);

/// Subtract a quaternion number from a real number
template <typename _F>
inline quaternion<_F> operator-(const _F&, const quaternion<_F>&);

/// Right multiply a quaternion number by a real number
template <typename _F>
inline quaternion<_F> operator*(const quaternion<_F>&, const _F&);

/// Left multiply a quaternion number by a real number
template <typename _F>
inline quaternion<_F> operator*(const _F&, const quaternion<_F>&);

/// Divide a quaternion number by a real number
template <typename _F>
inline quaternion<_F> operator/(const quaternion<_F>&, const _F&);

/// Divide a real number by a quaternion number
template <typename _F>
inline quaternion<_F> operator*(const _F&, const quaternion<_F>&);

/* @} */ // QROp

/**
 *  @defgroup QCOp Quaternion - Complex Binary Operations
 *  Binary arithmetic operations between quaternion and complex numbers
 *  @{
 */

/// Add a complex number to a quaternion number
template <typename _F>
inline quaternion<_F> operator+(const quaternion<_F>&, const std::complex<_F>&);

/// Add a complex number to a quaternion number
template <typename _F>
inline quaternion<_F> operator+(const std::complex<_F>&, const quaternion<_F>&);

/// Subtract a complex number from a quaternion number
template <typename _F>
inline quaternion<_F> operator-(const quaternion<_F>&, const std::complex<_F>&);

/// Subtract a quaternion number from a complex number
template <typename _F>
inline quaternion<_F> operator-(const std::complex<_F>&, const quaternion<_F>&);

/// Right multiply a quaternion number by a complex number
template <typename _F>
inline quaternion<_F> operator*(const quaternion<_F>&, const std::complex<_F>&);

/// Left multiply a quaternion number by a complex number
template <typename _F>
inline quaternion<_F> operator*(const std::complex<_F>&, const quaternion<_F>&);

/* @} */ // QCOp

/**
 *  @defgroup QAlg Quaternion Algebraic Operations
 *  Defines the algebraic operations for the quaternion numbers
 *  @{
 */

/// Returns the norm of a quaternion number
template <typename _F> inline _F norm(const quaternion<_F>&);

/// Returns the conjugate of a quaternion number
template <typename _F> inline quaternion<_F> conj(const quaternion<_F>&);

/// Returns the inverse of a quaternion number
template <typename _F> inline quaternion<_F>  inv(const quaternion<_F>&);

/// Returns the versor of a quaternion number
template <typename _F> inline quaternion<_F> versor(const quaternion<_F>&);

/// Returns the commutator of two quaternion numbers
template <typename _F> 
  inline quaternion<_F> comm(const quaternion<_F>&, const quaternion<_F>&);

/// Returns the commutator of a quaternion number and a complex number
template <typename _F> 
  inline quaternion<_F> comm(const quaternion<_F>&, const std::complex<_F>&);

/// Returns the commutator of a quaternion number and a complex number
template <typename _F> 
  inline quaternion<_F> comm(const std::complex<_F>&, const quaternion<_F>&);

/// Returns the commutator of a quaternion number and a real number
template <typename _F> 
  inline quaternion<_F> comm(const quaternion<_F>&, const _F&);

/// Returns the commutator of a quaternion number and a real number
template <typename _F> 
  inline quaternion<_F> comm(const _F&, const quaternion<_F>&);

/* @} */ // QAlg

/* @} */ // quaternion



/**
 *  @defgroup MiscFunctions Miscellaneous Function Overloads and Replacements
 *  Various overloads / replacements of standard functions
 *
 *  @{
 */

template <typename _F> inline _F SmartConj( _F& );
template<> inline double SmartConj( double &x ) { return x; }
template<> 
  inline std::complex<double> SmartConj( std::complex<double> &x ) { 
    return std::conj(x); 
  }
template<> 
  inline quaternion<double> SmartConj( quaternion<double> &x ) { 
    return conj(x); 
  }

/* @} */ // MiscFunctions 

}; // HAXX namespace



#endif
