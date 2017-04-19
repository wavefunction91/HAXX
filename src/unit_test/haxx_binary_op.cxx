#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HAXX_BINARY_OP

#ifndef _HAXX_UT_BUTF_NINCLUDED
  #include <boost/test/included/unit_test.hpp>
#else
  #include <boost/test/unit_test.hpp>
#endif

#include "haxx.hpp"


BOOST_AUTO_TEST_CASE(real_binary_add)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  double x = -4.59;
  p = q + x;
  r = x + q;

  BOOST_CHECK_EQUAL(p.real(),  q.real() + x);
  BOOST_CHECK_EQUAL(p.imag_i(),q.imag_i()  );
  BOOST_CHECK_EQUAL(p.imag_j(),q.imag_j()  );
  BOOST_CHECK_EQUAL(p.imag_k(),q.imag_k()  );

  BOOST_CHECK_EQUAL(r.real(),  q.real() + x);
  BOOST_CHECK_EQUAL(r.imag_i(),q.imag_i()  );
  BOOST_CHECK_EQUAL(r.imag_j(),q.imag_j()  );
  BOOST_CHECK_EQUAL(r.imag_k(),q.imag_k()  );

  BOOST_CHECK_EQUAL(q.real(),  1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);
}


BOOST_AUTO_TEST_CASE(real_binary_sub)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  double x = -4.59;
  p = q - x;
  r = x - q;

  BOOST_CHECK_EQUAL(p.real(),  q.real() - x);
  BOOST_CHECK_EQUAL(p.imag_i(),q.imag_i()  );
  BOOST_CHECK_EQUAL(p.imag_j(),q.imag_j()  );
  BOOST_CHECK_EQUAL(p.imag_k(),q.imag_k()  );

  BOOST_CHECK_EQUAL(r.real(),  x - q.real());
  BOOST_CHECK_EQUAL(r.imag_i(),-q.imag_i() );
  BOOST_CHECK_EQUAL(r.imag_j(),-q.imag_j() );
  BOOST_CHECK_EQUAL(r.imag_k(),-q.imag_k() );

  BOOST_CHECK_EQUAL(q.real(),  1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);
}

BOOST_AUTO_TEST_CASE(real_binary_mul)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  double x = -4.59;
  p = q * x;
  r = x * q;

  BOOST_CHECK_EQUAL(p.real(),  q.real() * x);
  BOOST_CHECK_EQUAL(p.imag_i(),q.imag_i()*x);
  BOOST_CHECK_EQUAL(p.imag_j(),q.imag_j()*x);
  BOOST_CHECK_EQUAL(p.imag_k(),q.imag_k()*x);

  BOOST_CHECK_EQUAL(r.real(),  q.real() * x);
  BOOST_CHECK_EQUAL(r.imag_i(),q.imag_i() *x);
  BOOST_CHECK_EQUAL(r.imag_j(),q.imag_j() *x);
  BOOST_CHECK_EQUAL(r.imag_k(),q.imag_k() *x);

  BOOST_CHECK_EQUAL(q.real(),  1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);
}

BOOST_AUTO_TEST_CASE(real_binary_div)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  double x = -4.59;
  p = q / x;
  r = x / q;

  double nrm = HAXX::norm(q);
  nrm *= nrm;

  BOOST_CHECK_EQUAL(p.real(),  q.real() / x);
  BOOST_CHECK_EQUAL(p.imag_i(),q.imag_i()/x);
  BOOST_CHECK_EQUAL(p.imag_j(),q.imag_j()/x);
  BOOST_CHECK_EQUAL(p.imag_k(),q.imag_k()/x);

  BOOST_CHECK_EQUAL(r.real(),   x*q.real()  /nrm);
  BOOST_CHECK_EQUAL(r.imag_i(),-x*q.imag_i()/nrm);
  BOOST_CHECK_EQUAL(r.imag_j(),-x*q.imag_j()/nrm);
  BOOST_CHECK_EQUAL(r.imag_k(),-x*q.imag_k()/nrm);

  BOOST_CHECK_EQUAL(q.real(),  1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);
}



BOOST_AUTO_TEST_CASE(complex_binary_add)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  std::complex<double> x(-4.59,4.928364);
  p = q + x;
  r = x + q;

  BOOST_CHECK_EQUAL(p.real(),  q.real() + x.real());
  BOOST_CHECK_EQUAL(p.imag_i(),q.imag_i() + x.imag()  );
  BOOST_CHECK_EQUAL(p.imag_j(),q.imag_j()  );
  BOOST_CHECK_EQUAL(p.imag_k(),q.imag_k()  );

  BOOST_CHECK_EQUAL(r.real(),  q.real() + x.real());
  BOOST_CHECK_EQUAL(r.imag_i(),q.imag_i() + x.imag() );
  BOOST_CHECK_EQUAL(r.imag_j(),q.imag_j()  );
  BOOST_CHECK_EQUAL(r.imag_k(),q.imag_k()  );

  BOOST_CHECK_EQUAL(q.real(),  1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);
}


BOOST_AUTO_TEST_CASE(complex_binary_sub)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.),p,r;
  std::complex<double> x(-4.59,4.928364);
  p = q - x;
  r = x - q;

  BOOST_CHECK_EQUAL(p.real(),  q.real() - x.real());
  BOOST_CHECK_EQUAL(p.imag_i(),q.imag_i()-x.imag()  );
  BOOST_CHECK_EQUAL(p.imag_j(),q.imag_j()  );
  BOOST_CHECK_EQUAL(p.imag_k(),q.imag_k()  );

  BOOST_CHECK_EQUAL(r.real(),  x.real() - q.real());
  BOOST_CHECK_EQUAL(r.imag_i(),x.imag() - q.imag_i() );
  BOOST_CHECK_EQUAL(r.imag_j(),-q.imag_j() );
  BOOST_CHECK_EQUAL(r.imag_k(),-q.imag_k() );

  BOOST_CHECK_EQUAL(q.real(),  1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);
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

  BOOST_CHECK_EQUAL(p.real()  ,S);
  BOOST_CHECK_EQUAL(p.imag_i(),I);
  BOOST_CHECK_EQUAL(p.imag_j(),J);
  BOOST_CHECK_EQUAL(p.imag_k(),K);

  S = B1 * A1 - B2 * A2;
  I = B1 * A2 + B2 * A1;
  J = B1 * A3 - B2 * A4;
  K = B1 * A4 + B2 * A3;

  BOOST_CHECK_EQUAL(r.real()  ,S);
  BOOST_CHECK_EQUAL(r.imag_i(),I);
  BOOST_CHECK_EQUAL(r.imag_j(),J);
  BOOST_CHECK_EQUAL(r.imag_k(),K);


  BOOST_CHECK_EQUAL(q.real(),  1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);
}


BOOST_AUTO_TEST_CASE(quaternion_binary_add)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.), p(-5.,6.94,7.82,-8.),r;
  r = p + q;

  BOOST_CHECK_EQUAL(r.real()  ,p.real()   + q.real()  );
  BOOST_CHECK_EQUAL(r.imag_i(),p.imag_i() + q.imag_i());
  BOOST_CHECK_EQUAL(r.imag_j(),p.imag_j() + q.imag_j());
  BOOST_CHECK_EQUAL(r.imag_k(),p.imag_k() + q.imag_k());

  BOOST_CHECK_EQUAL(q.real(),1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);

  BOOST_CHECK_EQUAL(p.real()  ,-5.);
  BOOST_CHECK_EQUAL(p.imag_i(),6.94);
  BOOST_CHECK_EQUAL(p.imag_j(),7.82);
  BOOST_CHECK_EQUAL(p.imag_k(),-8.);


}


BOOST_AUTO_TEST_CASE(quaternion_binary_sub)
{

  HAXX::quaternion<double> q(1.,2.,3.,4.), p(-5.,6.94,7.82,-8.),r;
  r = p - q;

  BOOST_CHECK_EQUAL(r.real()  ,p.real()   - q.real()  );
  BOOST_CHECK_EQUAL(r.imag_i(),p.imag_i() - q.imag_i());
  BOOST_CHECK_EQUAL(r.imag_j(),p.imag_j() - q.imag_j());
  BOOST_CHECK_EQUAL(r.imag_k(),p.imag_k() - q.imag_k());

  BOOST_CHECK_EQUAL(q.real(),1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);

  BOOST_CHECK_EQUAL(p.real()  ,-5.);
  BOOST_CHECK_EQUAL(p.imag_i(),6.94);
  BOOST_CHECK_EQUAL(p.imag_j(),7.82);
  BOOST_CHECK_EQUAL(p.imag_k(),-8.);


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

  BOOST_CHECK_EQUAL(r.real()  ,S);
  BOOST_CHECK_EQUAL(r.imag_i(),I);
  BOOST_CHECK_EQUAL(r.imag_j(),J);
  BOOST_CHECK_EQUAL(r.imag_k(),K);

  BOOST_CHECK_EQUAL(q.real(),1.);
  BOOST_CHECK_EQUAL(q.imag_i(),2.);
  BOOST_CHECK_EQUAL(q.imag_j(),3.);
  BOOST_CHECK_EQUAL(q.imag_k(),4.);

  BOOST_CHECK_EQUAL(p.real()  ,-5.);
  BOOST_CHECK_EQUAL(p.imag_i(),6.94);
  BOOST_CHECK_EQUAL(p.imag_j(),7.82);
  BOOST_CHECK_EQUAL(p.imag_k(),-8.);


}



