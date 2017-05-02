/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HAXX_ALGEBRA

#include "haxx_ut.hpp"


BOOST_AUTO_TEST_CASE(conj)
{
  HAXX::quaternion<double> q(1.,2.,3.,4.);
  HAXX::quaternion<double> p = HAXX::conj(q);

  
  BOOST_CHECK_EQUAL(p.real(),  q.real()   );
  BOOST_CHECK_EQUAL(p.imag_i(),-q.imag_i());
  BOOST_CHECK_EQUAL(p.imag_j(),-q.imag_j());
  BOOST_CHECK_EQUAL(p.imag_k(),-q.imag_k());

  BOOST_CHECK_EQUAL(q.real(),  1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);
};

BOOST_AUTO_TEST_CASE(norm)
{
  HAXX::quaternion<double> q(1.,2.,3.,4.);
  double nrm = HAXX::norm(q);
  double nrmStupid = 
    q.real() * q.real() +
    q.imag_i() * q.imag_i() +
    q.imag_j() * q.imag_j() +
    q.imag_k() * q.imag_k();

  nrmStupid = std::sqrt(nrmStupid);


  HAXX::quaternion<double> p(8.43,-9.3829,-10.3,3.14);
  HAXX::quaternion<double> q1(7.,2.43,3.65,7.006);
  HAXX::quaternion<double> p1(4.,2.,3.,7.);
  double a = 0.543;

  // Distance function
  double distPQ = HAXX::norm(p - q);
  double distQP = HAXX::norm(q - p);

  // Scaled distance
  double distP1Q1 = HAXX::norm(p1 + q1);
  double distDiff = HAXX::norm((p + a*p1 + q + a*q1) - (p + q));

  
  BOOST_CHECK_EQUAL(nrm,nrmStupid);
  BOOST_CHECK_EQUAL(distPQ,distQP);
  BOOST_CHECK_EQUAL(distDiff,a*distP1Q1);

  BOOST_CHECK_EQUAL(q.real(),  1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);

};

BOOST_AUTO_TEST_CASE(hpow)
{
  double a = genRandom<double>();
  HAXX::quaternion<double> q = genRandom<HAXX::quaternion<double>>();

  HAXX::quaternion<double> qa = HAXX::pow(q,a);
  HAXX::quaternion<double> qma = HAXX::pow(q,-a);

  HAXX::quaternion<double> prod1 = qa * qma;
  HAXX::quaternion<double> prod2 = qma * qa;

  BOOST_CHECK( CMP_Q(prod1,prod2) );
  BOOST_CHECK( CMP_Q(prod1,HAXX::quaternion<double>(1.))    );
};

BOOST_AUTO_TEST_CASE(hsqrt)
{
  HAXX::quaternion<double> q = genRandom<HAXX::quaternion<double>>();
  HAXX::quaternion<double> qr = HAXX::sqrt(q);
  HAXX::quaternion<double> prod1 = qr * qr;

  BOOST_CHECK( CMP_Q(prod1,q) );
};
