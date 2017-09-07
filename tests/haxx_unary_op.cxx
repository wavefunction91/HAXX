/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */


#include "haxx_ut.hpp"
 
BOOST_AUTO_TEST_SUITE(HAXX_UNARY_OP)

// Real Unary Operators
BOOST_AUTO_TEST_CASE(real_unary_assign)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.);
  double x = -4.59;
  q = x;

  // Check that q is unchanged
  BOOST_CHECK_EQUAL(q.real(),x);
  BOOST_CHECK_EQUAL(q.imag_i(),0.);
  BOOST_CHECK_EQUAL(q.imag_j(),0.);
  BOOST_CHECK_EQUAL(q.imag_k(),0.);

}

BOOST_AUTO_TEST_CASE(real_unary_add)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.);
  double x = -4.59;
  q += x;

  // Check that q is unchanged
  BOOST_CHECK_EQUAL(q.real(),1. + (x));
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);

}


BOOST_AUTO_TEST_CASE(real_unary_sub)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.);
  double x = -4.59;
  q -= x;

  // Check that q is unchanged
  BOOST_CHECK_EQUAL(q.real(),1. - (x));
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);

}

BOOST_AUTO_TEST_CASE(real_unary_mul)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.);
  double x = -4.59;
  q *= x;

  // Check that q is unchanged
  BOOST_CHECK_EQUAL(q.real(),  1.*x);
  BOOST_CHECK_EQUAL(q.imag_i(),2.*x);
  BOOST_CHECK_EQUAL(q.imag_j(),3.*x);
  BOOST_CHECK_EQUAL(q.imag_k(),4.*x);

}

BOOST_AUTO_TEST_CASE(real_unary_div)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.);
  double x = -4.59;
  q /= x;

  // Check that q is unchanged
  BOOST_CHECK_EQUAL(q.real(),  1./x);
  BOOST_CHECK_EQUAL(q.imag_i(),2./x);
  BOOST_CHECK_EQUAL(q.imag_j(),3./x);
  BOOST_CHECK_EQUAL(q.imag_k(),4./x);

}



// Complex Unary Operators
BOOST_AUTO_TEST_CASE(complex_unary_assign)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.);
  std::complex<double> x(-4.76,5.6);
  q = x;

  // Check that q is unchanged
  BOOST_CHECK_EQUAL(q.real(),x.real());
  BOOST_CHECK_EQUAL(q.imag_i(),x.imag());
  BOOST_CHECK_EQUAL(q.imag_j(),0.);
  BOOST_CHECK_EQUAL(q.imag_k(),0.);

}

BOOST_AUTO_TEST_CASE(complex_unary_add)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.);
  std::complex<double> x(-4.76,5.6);
  q += x;

  // Check that q is unchanged
  BOOST_CHECK_EQUAL(q.real(),1. + x.real());
  BOOST_CHECK_EQUAL(q.imag_i(),2. + x.imag());
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);

}


BOOST_AUTO_TEST_CASE(complex_unary_sub)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.);
  std::complex<double> x(-4.76,5.6);
  q -= x;

  // Check that q is unchanged
  BOOST_CHECK_EQUAL(q.real(),1. - x.real());
  BOOST_CHECK_EQUAL(q.imag_i(),2. - x.imag());
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);

}


// Quaternion Unary Operations
BOOST_AUTO_TEST_CASE(quaternion_unary_assign)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.), p = q;

  // Check that q is unchanged
  BOOST_CHECK_EQUAL(q.real(),p.real());
  BOOST_CHECK_EQUAL(q.imag_i(),p.imag_i());
  BOOST_CHECK_EQUAL(q.imag_j(),p.imag_j());
  BOOST_CHECK_EQUAL(q.imag_k(),p.imag_k());

}

BOOST_AUTO_TEST_CASE(quaternion_unary_add)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.), p(5.,6.,7.,8.);
  p += q;

  // Check that p has been changed properly
  BOOST_CHECK_EQUAL(p.real(),6.);
  BOOST_CHECK_EQUAL(p.imag_i(),8.);
  BOOST_CHECK_EQUAL(p.imag_j(),10.);
  BOOST_CHECK_EQUAL(p.imag_k(),12.);

  // Check that q is unchanged
  BOOST_CHECK_EQUAL(q.real(),1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);
}

BOOST_AUTO_TEST_CASE(quaternion_unary_sub)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.), p(5.,0.,7.,9.);
  p -= q;

  // Check that p has been changed properly
  BOOST_CHECK_EQUAL(p.real(),4.);
  BOOST_CHECK_EQUAL(p.imag_i(),-2.);
  BOOST_CHECK_EQUAL(p.imag_j(),4.);
  BOOST_CHECK_EQUAL(p.imag_k(),5.);

  // Check that q is unchanged
  BOOST_CHECK_EQUAL(q.real(),1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);

}

BOOST_AUTO_TEST_SUITE_END()
