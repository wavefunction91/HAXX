#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HLAPACK_HOUSEHOLDER

#include "haxx_ut.hpp"
#include "hblas/hblas1_impl.hpp"
#include "hlapack/householder_impl.hpp"


BOOST_AUTO_TEST_CASE(hlapack_larfg)
{
  std::vector<HAXX::quaternion<double>> 
    X(HBLAS1_VECLEN), XC(HBLAS1_VECLEN), SCR(HBLAS1_VECLEN);

  
  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

/*
  for(auto &x : X) 
    std::cout << x << std::endl;
  std::cout << std::endl; 
*/

  XC = X;

  double BETA = HLAPACK_LARFG(HBLAS1_VECLEN,&X[0],1);
  
  auto Inner = HBLAS_DOTC(HBLAS1_VECLEN,&X[0],1,&XC[0],1);
  HBLAS_SCAL('R',HBLAS1_VECLEN,Inner,&X[0],1);
  HBLAS_AXPY('L',HBLAS1_VECLEN,-BETA,&X[0],1,&XC[0],1);

/*
  for(auto &x : XC) 
    std::cout << x << std::endl;
  std::cout << std::endl; 
*/

  double x2nrm = HBLAS_NRM2(HBLAS1_VECLEN,&XC[0],1);
  auto   Uy    = XC[0] / HAXX::norm(XC[0]);

  BOOST_CHECK(CMP_Q(x2nrm*Uy,XC[0]));
  for(auto j = 1; j < HBLAS1_VECLEN; ++j){
    BOOST_CHECK_SMALL(XC[j].real()  ,COMPARE_TOL);
    BOOST_CHECK_SMALL(XC[j].imag_i(),COMPARE_TOL);
    BOOST_CHECK_SMALL(XC[j].imag_j(),COMPARE_TOL);
    BOOST_CHECK_SMALL(XC[j].imag_k(),COMPARE_TOL);
  }

}
