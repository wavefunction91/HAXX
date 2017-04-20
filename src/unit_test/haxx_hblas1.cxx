#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HAXX_BLAS1

#ifndef _HAXX_UT_BUTF_NINCLUDED
  #include <boost/test/included/unit_test.hpp>
#else
  #include <boost/test/unit_test.hpp>
#endif

#include "haxx.hpp"
#include "hblas/hblas1_def.hpp"

#include <random>
#include <iterator>
#include <iostream>

#define HBLAS1_VECLEN 500
#define HBLAS1_RAND_MIN -20
#define HBLAS1_RAND_MAX 54

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(HBLAS1_RAND_MIN,HBLAS1_RAND_MAX);

BOOST_AUTO_TEST_CASE(hblas1_swap)
{
  // Random Quaternion vector
  std::vector<HAXX::quaternion<double>> X(HBLAS1_VECLEN), Y(HBLAS1_VECLEN);
  for(auto &x : X) x = dis(gen);
  for(auto &x : Y) x = dis(gen);

  std::vector<HAXX::quaternion<double>> XC(HBLAS1_VECLEN),YC(HBLAS1_VECLEN);
  std::copy(X.begin(),X.end(),XC.begin());
  std::copy(Y.begin(),Y.end(),YC.begin());

  SWAP(HBLAS1_VECLEN,&X[0],1,&Y[0],1);
}

BOOST_AUTO_TEST_CASE(hblas1_scal)
{

  // Random Quaternion vector
  std::vector<HAXX::quaternion<double>> X(HBLAS1_VECLEN);
  for(auto &x : X) x = dis(gen);

  // Random scaling factors
  double                   rAlpha = dis(gen);
  std::complex<double>     cAlpha(dis(gen),dis(gen));
  HAXX::quaternion<double> hAlpha(dis(gen),dis(gen),dis(gen),dis(gen));
 
  std::cout << "hblas_scal will test the following scaling factors:" <<
    std::endl;
  std::cout << "  real : " << rAlpha << std::endl;
  std::cout << "  complex : " << cAlpha << std::endl;
  std::cout << "  quaternion : " << hAlpha << std::endl;

  std::vector<HAXX::quaternion<double>> tmpX(HBLAS1_VECLEN);
  std::vector<size_t> strides = {1,2,3,5,9};

  // Right scaling by real alpha
  for(auto stride : strides) {
    std::copy(tmpX.begin(),tmpX.end(),X.begin());
    SCAL('R',HBLAS1_VECLEN,rAlpha,&tmpX[0],stride);
    BOOST_CHECK( 
      std::all_of(tmpX.begin(),tmpX.end(),
        [&](HAXX::quaternion<double> &x) {
          size_t indx = std::distance(&tmpX[0],&x);
	  if(indx % stride == 0) return (x == X[indx]*rAlpha);
	  else return (x == X[indx]);
        }
      )
    );
  }

  // Left scaling by real alpha
  for(auto stride : strides) {
    std::copy(tmpX.begin(),tmpX.end(),X.begin());
    SCAL('L',HBLAS1_VECLEN,rAlpha,&tmpX[0],stride);
    BOOST_CHECK( 
      std::all_of(tmpX.begin(),tmpX.end(),
        [&](HAXX::quaternion<double> &x) {
          size_t indx = std::distance(&tmpX[0],&x);
	  if(indx % stride == 0) return (x == X[indx]*rAlpha);
	  else return (x == X[indx]);
        }
      )
    );
  }

  // Right scaling by complex alpha
  for(auto stride : strides) {
    std::copy(tmpX.begin(),tmpX.end(),X.begin());
    SCAL('R',HBLAS1_VECLEN,cAlpha,&tmpX[0],stride);
    BOOST_CHECK( 
      std::all_of(tmpX.begin(),tmpX.end(),
        [&](HAXX::quaternion<double> &x) {
          size_t indx = std::distance(&tmpX[0],&x);
	  if(indx % stride == 0) return (x == X[indx]*cAlpha);
	  else return (x == X[indx]);
        }
      )
    );
  }

  // Left scaling by complex alpha
  for(auto stride : strides) {
    std::copy(tmpX.begin(),tmpX.end(),X.begin());
    SCAL('L',HBLAS1_VECLEN,cAlpha,&tmpX[0],stride);
    BOOST_CHECK( 
      std::all_of(tmpX.begin(),tmpX.end(),
        [&](HAXX::quaternion<double> &x) {
          size_t indx = std::distance(&tmpX[0],&x);
	  if(indx % stride == 0) return (x == cAlpha*X[indx]);
	  else return (x == X[indx]);
        }
      )
    );
  }

  // Right scaling by quaternion alpha
  for(auto stride : strides) {
    std::copy(tmpX.begin(),tmpX.end(),X.begin());
    SCAL('R',HBLAS1_VECLEN,hAlpha,&tmpX[0],stride);
    BOOST_CHECK( 
      std::all_of(tmpX.begin(),tmpX.end(),
        [&](HAXX::quaternion<double> &x) {
          size_t indx = std::distance(&tmpX[0],&x);
	  if(indx % stride == 0) return (x == X[indx]*hAlpha);
	  else return (x == X[indx]);
        }
      )
    );
  }

  // Left scaling by quaternion alpha
  for(auto stride : strides) {
    std::copy(tmpX.begin(),tmpX.end(),X.begin());
    SCAL('L',HBLAS1_VECLEN,hAlpha,&tmpX[0],stride);
    BOOST_CHECK( 
      std::all_of(tmpX.begin(),tmpX.end(),
        [&](HAXX::quaternion<double> &x) {
          size_t indx = std::distance(&tmpX[0],&x);
	  if(indx % stride == 0) return (x == hAlpha*X[indx]);
	  else return (x == X[indx]);
        }
      )
    );
  }
};
