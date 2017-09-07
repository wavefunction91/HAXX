/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */


#include "haxx_ut.hpp"

BOOST_AUTO_TEST_SUITE(HAXX_BINARY_OP)

BOOST_AUTO_TEST_CASE(real_binary_add)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  double x = -4.59;
  p = q + x;
  r = x + q;

  BOOST_CHECK_CLOSE(p.real(),  q.real() + x, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_i(),q.imag_i()  , COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_j(),q.imag_j()  , COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_k(),q.imag_k()  , COMPARE_TOL);

  BOOST_CHECK_CLOSE(r.real(),  q.real() + x, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_i(),q.imag_i()  , COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_j(),q.imag_j()  , COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_k(),q.imag_k()  , COMPARE_TOL);

  BOOST_CHECK_CLOSE(q.real(),  1., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_i(),2., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_j(),3., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_k(),4., COMPARE_TOL);
}


BOOST_AUTO_TEST_CASE(real_binary_sub)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  double x = -4.59;
  p = q - x;
  r = x - q;

  BOOST_CHECK_CLOSE(p.real(),  q.real() - x, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_i(),q.imag_i()  , COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_j(),q.imag_j()  , COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_k(),q.imag_k()  , COMPARE_TOL);

  BOOST_CHECK_CLOSE(r.real(),  x - q.real(), COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_i(),-q.imag_i() , COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_j(),-q.imag_j() , COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_k(),-q.imag_k() , COMPARE_TOL);

  BOOST_CHECK_CLOSE(q.real(),  1., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_i(),2., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_j(),3., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_k(),4., COMPARE_TOL);
}

BOOST_AUTO_TEST_CASE(real_binary_mul)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  double x = -4.59;
  p = q * x;
  r = x * q;

  BOOST_CHECK_CLOSE(p.real(),  q.real() * x, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_i(),q.imag_i()*x, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_j(),q.imag_j()*x, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_k(),q.imag_k()*x, COMPARE_TOL);

  BOOST_CHECK_CLOSE(r.real(),  q.real() * x, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_i(),q.imag_i() *x, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_j(),q.imag_j() *x, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_k(),q.imag_k() *x, COMPARE_TOL);

  BOOST_CHECK_CLOSE(q.real(),  1., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_i(),2., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_j(),3., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_k(),4., COMPARE_TOL);
}

BOOST_AUTO_TEST_CASE(real_binary_div)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  double x = -4.59;
  p = q / x;
  r = x / q;

  double nrm = HAXX::norm(q);
  nrm *= nrm;

  BOOST_CHECK_CLOSE(p.real(),  q.real() / x, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_i(),q.imag_i()/x, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_j(),q.imag_j()/x, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_k(),q.imag_k()/x, COMPARE_TOL);

  BOOST_CHECK_CLOSE(r.real(),   x*q.real()  /nrm, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_i(),-x*q.imag_i()/nrm, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_j(),-x*q.imag_j()/nrm, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_k(),-x*q.imag_k()/nrm, COMPARE_TOL);

  BOOST_CHECK_CLOSE(q.real(),  1., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_i(),2., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_j(),3., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_k(),4., COMPARE_TOL);
}



BOOST_AUTO_TEST_CASE(complex_binary_add)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  std::complex<double> x(-4.59,4.928364);
  p = q + x;
  r = x + q;

  BOOST_CHECK_CLOSE(p.real(),  q.real() + x.real()  , COMPARE_TOL  );
  BOOST_CHECK_CLOSE(p.imag_i(),q.imag_i() + x.imag(), COMPARE_TOL  );
  BOOST_CHECK_CLOSE(p.imag_j(),q.imag_j()           , COMPARE_TOL  );
  BOOST_CHECK_CLOSE(p.imag_k(),q.imag_k()           , COMPARE_TOL  );

  BOOST_CHECK_CLOSE(r.real(),  q.real() + x.real()  , COMPARE_TOL );
  BOOST_CHECK_CLOSE(r.imag_i(),q.imag_i() + x.imag(), COMPARE_TOL );
  BOOST_CHECK_CLOSE(r.imag_j(),q.imag_j()           , COMPARE_TOL );
  BOOST_CHECK_CLOSE(r.imag_k(),q.imag_k()           , COMPARE_TOL );

  BOOST_CHECK_CLOSE(q.real(),  1., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_i(),2., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_j(),3., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_k(),4., COMPARE_TOL);
}


BOOST_AUTO_TEST_CASE(complex_binary_sub)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  std::complex<double> x(-4.59,4.928364);
  p = q - x;
  r = x - q;

  BOOST_CHECK_CLOSE(p.real(),  q.real() - x.real(), COMPARE_TOL  );
  BOOST_CHECK_CLOSE(p.imag_i(),q.imag_i()-x.imag(), COMPARE_TOL  );
  BOOST_CHECK_CLOSE(p.imag_j(),q.imag_j()         , COMPARE_TOL  );
  BOOST_CHECK_CLOSE(p.imag_k(),q.imag_k()         , COMPARE_TOL  );

  BOOST_CHECK_CLOSE(r.real(),  x.real() - q.real()  , COMPARE_TOL );
  BOOST_CHECK_CLOSE(r.imag_i(),x.imag() - q.imag_i(), COMPARE_TOL );
  BOOST_CHECK_CLOSE(r.imag_j(),-q.imag_j()          , COMPARE_TOL );
  BOOST_CHECK_CLOSE(r.imag_k(),-q.imag_k()          , COMPARE_TOL );

  BOOST_CHECK_CLOSE(q.real(),  1., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_i(),2., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_j(),3., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_k(),4., COMPARE_TOL);
}



BOOST_AUTO_TEST_CASE(complex_binary_mul)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  std::complex<double> x(-4.59,4.928364);
  p = q * x;
  r = x * q;

  double A1 = q.real(), A2 = q.imag_i(), A3 = q.imag_j(), A4 = q.imag_k();
  double B1 = x.real(), B2 = x.imag();

  double S = A1 * B1 - A2 * B2;
  double I = A1 * B2 + A2 * B1;
  double J = A3 * B1 + A4 * B2;
  double K = A4 * B1 - A3 * B2;

  BOOST_CHECK_CLOSE(p.real()  ,S, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_i(),I, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_j(),J, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_k(),K, COMPARE_TOL);

  S = B1 * A1 - B2 * A2;
  I = B1 * A2 + B2 * A1;
  J = B1 * A3 - B2 * A4;
  K = B1 * A4 + B2 * A3;

  BOOST_CHECK_CLOSE(r.real()  ,S, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_i(),I, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_j(),J, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_k(),K, COMPARE_TOL);


  BOOST_CHECK_CLOSE(q.real(),  1., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_i(),2., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_j(),3., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_k(),4., COMPARE_TOL);
}


BOOST_AUTO_TEST_CASE(quaternion_binary_add)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.), p(-5.,6.94,7.82,-8.),r;
  r = p + q;

  BOOST_CHECK_CLOSE(r.real()  ,p.real()   + q.real()  , COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_i(),p.imag_i() + q.imag_i(), COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_j(),p.imag_j() + q.imag_j(), COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_k(),p.imag_k() + q.imag_k(), COMPARE_TOL);

  BOOST_CHECK_CLOSE(q.real(),1.  , COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_i(),2., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_j(),3., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_k(),4., COMPARE_TOL);

  BOOST_CHECK_CLOSE(p.real()  ,-5. , COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_i(),6.94, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_j(),7.82, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_k(),-8. , COMPARE_TOL);


}


BOOST_AUTO_TEST_CASE(quaternion_binary_sub)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.), p(-5.,6.94,7.82,-8.),r;
  r = p - q;

  BOOST_CHECK_CLOSE(r.real()  ,p.real()   - q.real()  , COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_i(),p.imag_i() - q.imag_i(), COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_j(),p.imag_j() - q.imag_j(), COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_k(),p.imag_k() - q.imag_k(), COMPARE_TOL);

  BOOST_CHECK_CLOSE(q.real(),1.  , COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_i(),2., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_j(),3., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_k(),4., COMPARE_TOL);

  BOOST_CHECK_CLOSE(p.real()  ,-5. , COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_i(),6.94, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_j(),7.82, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_k(),-8. , COMPARE_TOL);


}


BOOST_AUTO_TEST_CASE(quaternion_binary_mul)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.), p(-5.,6.94,7.82,-8.),r;
  r = p * q;

  double A1 = p.real(), A2 = p.imag_i(), A3 = p.imag_j(), A4 = p.imag_k();
  double B1 = q.real(), B2 = q.imag_i(), B3 = q.imag_j(), B4 = q.imag_k();

  double S = A1 * B1 - A2 * B2 - A3 * B3 - A4 * B4;
  double I = A1 * B2 + A2 * B1 + A3 * B4 - A4 * B3;
  double J = A1 * B3 - A2 * B4 + A3 * B1 + A4 * B2;
  double K = A1 * B4 + A2 * B3 - A3 * B2 + A4 * B1;

  BOOST_CHECK_CLOSE(r.real()  ,S, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_i(),I, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_j(),J, COMPARE_TOL);
  BOOST_CHECK_CLOSE(r.imag_k(),K, COMPARE_TOL);

  BOOST_CHECK_CLOSE(q.real(),1.  , COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_i(),2., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_j(),3., COMPARE_TOL);
  BOOST_CHECK_CLOSE(q.imag_k(),4., COMPARE_TOL);

  BOOST_CHECK_CLOSE(p.real()  ,-5. , COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_i(),6.94, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_j(),7.82, COMPARE_TOL);
  BOOST_CHECK_CLOSE(p.imag_k(),-8. , COMPARE_TOL);


}


BOOST_AUTO_TEST_SUITE_END()
