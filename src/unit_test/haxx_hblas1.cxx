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
  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : Y)
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  // Make copies
  std::vector<HAXX::quaternion<double>> XC(X),YC(Y);
  
  // Index list
  std::vector<int> indx(HBLAS1_VECLEN); 
  std::iota(indx.begin(),indx.end(),0);


  std::vector<size_t> strides = {1,2,3,5,9};

  // Check swap when both strides are equal
  for(int stride : strides) { 

    std::copy(XC.begin(),XC.end(),X.begin());
    std::copy(YC.begin(),YC.end(),Y.begin());

    int len = HBLAS1_VECLEN/stride;
    SWAP(len,&X[0],stride,&Y[0],stride);
  
    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
          if(indx % stride == 0 and indx != len*stride ) 
            return (X[indx] == YC[indx] and Y[indx] == XC[indx]);
          else
            return (X[indx] == XC[indx] and Y[indx] == YC[indx]);
        }
      )
    );

  }

  // FIXME: Need a test for when strides are not equal
}

BOOST_AUTO_TEST_CASE(hblas1_scal)
{

  // Random Quaternion vector
  std::vector<HAXX::quaternion<double>> X(HBLAS1_VECLEN);
  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

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
    std::copy(X.begin(),X.end(),tmpX.begin());

    auto len = HBLAS1_VECLEN/stride; 
    SCAL('R',len,rAlpha,&tmpX[0],stride);

    BOOST_CHECK( 
      std::all_of(tmpX.begin(),tmpX.end(),
        [&](HAXX::quaternion<double> &x) {
          size_t indx = std::distance(&tmpX[0],&x);
	  if(indx % stride == 0 and indx != len*stride ) 
            return (x == X[indx]*rAlpha);
	  else return (x == X[indx]);
        }
      )
    );

  }


  // Left scaling by real alpha
  for(auto stride : strides) {
    std::copy(X.begin(),X.end(),tmpX.begin());

    auto len = HBLAS1_VECLEN/stride; 
    SCAL('L',len,rAlpha,&tmpX[0],stride);

    BOOST_CHECK( 
      std::all_of(tmpX.begin(),tmpX.end(),
        [&](HAXX::quaternion<double> &x) {
          size_t indx = std::distance(&tmpX[0],&x);
	  if(indx % stride == 0 and indx != len*stride ) 
            return (x == rAlpha*X[indx]);
	  else return (x == X[indx]);
        }
      )
    );

  }

  // Right scaling by complex alpha
  for(auto stride : strides) {
    std::copy(X.begin(),X.end(),tmpX.begin());

    auto len = HBLAS1_VECLEN/stride; 
    SCAL('R',len,cAlpha,&tmpX[0],stride);

    BOOST_CHECK( 
      std::all_of(tmpX.begin(),tmpX.end(),
        [&](HAXX::quaternion<double> &x) {
          size_t indx = std::distance(&tmpX[0],&x);
	  if(indx % stride == 0 and indx != len*stride ) 
            return (x == X[indx]*cAlpha);
	  else return (x == X[indx]);
        }
      )
    );

  }


  // Left scaling by complex alpha
  for(auto stride : strides) {
    std::copy(X.begin(),X.end(),tmpX.begin());

    auto len = HBLAS1_VECLEN/stride; 
    SCAL('L',len,cAlpha,&tmpX[0],stride);

    BOOST_CHECK( 
      std::all_of(tmpX.begin(),tmpX.end(),
        [&](HAXX::quaternion<double> &x) {
          size_t indx = std::distance(&tmpX[0],&x);
	  if(indx % stride == 0 and indx != len*stride ) 
            return (x == cAlpha*X[indx]);
	  else return (x == X[indx]);
        }
      )
    );

  }

  // Right scaling by quaternion alpha
  for(auto stride : strides) {
    std::copy(X.begin(),X.end(),tmpX.begin());

    auto len = HBLAS1_VECLEN/stride; 
    SCAL('R',len,hAlpha,&tmpX[0],stride);

    BOOST_CHECK( 
      std::all_of(tmpX.begin(),tmpX.end(),
        [&](HAXX::quaternion<double> &x) {
          size_t indx = std::distance(&tmpX[0],&x);
	  if(indx % stride == 0 and indx != len*stride ) 
            return (x == X[indx]*hAlpha);
	  else return (x == X[indx]);
        }
      )
    );

  }


  // Left scaling by quaternion alpha
  for(auto stride : strides) {
    std::copy(X.begin(),X.end(),tmpX.begin());

    auto len = HBLAS1_VECLEN/stride; 
    SCAL('L',len,hAlpha,&tmpX[0],stride);

    BOOST_CHECK( 
      std::all_of(tmpX.begin(),tmpX.end(),
        [&](HAXX::quaternion<double> &x) {
          size_t indx = std::distance(&tmpX[0],&x);
	  if(indx % stride == 0 and indx != len*stride ) 
            return (x == hAlpha*X[indx]);
	  else return (x == X[indx]);
        }
      )
    );

  }

};

BOOST_AUTO_TEST_CASE(hblas1_copy)
{
  // Random Quaternion vector
  std::vector<HAXX::quaternion<double>> X(HBLAS1_VECLEN), Y(HBLAS1_VECLEN);
  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));


  // Index list
  std::vector<int> indx(HBLAS1_VECLEN); 
  std::iota(indx.begin(),indx.end(),0);



  std::vector<size_t> strides = {1,2,3,5,9};

  // Check swap when both strides are equal
  for(int stride : strides) { 

    std::fill(Y.begin(),Y.end(),HAXX::quaternion<double>(0));

    int len = HBLAS1_VECLEN/stride;
    COPY(len,&X[0],stride,&Y[0],stride);
  
    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
          if(indx % stride == 0 and indx != len*stride ) 
            return (X[indx] == Y[indx]);
          else
            return (Y[indx] == HAXX::quaternion<double>(0));
        }
      )
    );

  }

  // FIXME: Need a test for when strides are not equal
}

BOOST_AUTO_TEST_CASE(hblas1_axpy)
{
  // Random Quaternion vectors
  std::vector<HAXX::quaternion<double>> X(HBLAS1_VECLEN), Y(HBLAS1_VECLEN);
  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : Y)
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  // Random scaling factors
  double                   rAlpha = dis(gen);
  std::complex<double>     cAlpha(dis(gen),dis(gen));
  HAXX::quaternion<double> hAlpha(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas1_axpy will test the following scaling factors:" <<
    std::endl;
  std::cout << "  real : " << rAlpha << std::endl;
  std::cout << "  complex : " << cAlpha << std::endl;
  std::cout << "  quaternion : " << hAlpha << std::endl;

  // Index list
  std::vector<int> indx(HBLAS1_VECLEN); 
  std::iota(indx.begin(),indx.end(),0);


  std::vector<HAXX::quaternion<double>> YC(Y);
  std::vector<size_t> strides = {1,2,3,5,9};
  //std::vector<size_t> strides = {1};
  
  // Both strides equal
  for(auto stride : strides) {
    auto len = HBLAS1_VECLEN/stride; 

    // Right scaling by real alpha
    std::copy(YC.begin(),YC.end(),Y.begin());

    AXPY('R',len,rAlpha,&X[0],stride,&Y[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (Y[indx] == YC[indx] + X[indx]*rAlpha);
	  else return (Y[indx] == YC[indx]);
        }
      )
    );

    // Left scaling by real alpha
    std::copy(YC.begin(),YC.end(),Y.begin());

    AXPY('L',len,rAlpha,&X[0],stride,&Y[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (Y[indx] == YC[indx] + rAlpha*X[indx]);
	  else return (Y[indx] == YC[indx]);
        }
      )
    );

    // Right scaling by complex alpha
    std::copy(YC.begin(),YC.end(),Y.begin());

    AXPY('R',len,cAlpha,&X[0],stride,&Y[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (Y[indx] == YC[indx] + X[indx]*cAlpha);
	  else return (Y[indx] == YC[indx]);
        }
      )
    );

    // Left scaling by complex alpha
    std::copy(YC.begin(),YC.end(),Y.begin());

    AXPY('L',len,cAlpha,&X[0],stride,&Y[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (Y[indx] == YC[indx] + cAlpha*X[indx]);
	  else return (Y[indx] == YC[indx]);
        }
      )
    );


    // Right scaling by quaternion alpha
    std::copy(YC.begin(),YC.end(),Y.begin());

    AXPY('R',len,hAlpha,&X[0],stride,&Y[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (Y[indx] == YC[indx] + X[indx]*hAlpha);
	  else return (Y[indx] == YC[indx]);
        }
      )
    );

    // Left scaling by quaternion alpha
    std::copy(YC.begin(),YC.end(),Y.begin());

    AXPY('L',len,hAlpha,&X[0],stride,&Y[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (Y[indx] == YC[indx] + hAlpha*X[indx]);
	  else return (Y[indx] == YC[indx]);
        }
      )
    );
  }


  // FIXME: Need a test for when strides are not equal

};

BOOST_AUTO_TEST_CASE(hblas1_dot)
{
  // Random Quaternion vectors
  std::vector<HAXX::quaternion<double>> X(HBLAS1_VECLEN), Y(HBLAS1_VECLEN);
  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : Y)
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  // Index list
  std::vector<int> indx(HBLAS1_VECLEN); 
  std::iota(indx.begin(),indx.end(),0);

  std::vector<size_t> strides = {1,2,3,5,9};
  //std::vector<size_t> strides = {1};

  // Both strides equal
  for(int stride : strides) {
    auto len = HBLAS1_VECLEN/stride; 
    HAXX::quaternion<double> dotu = DOTU(len,&X[0],stride,&Y[0],stride);
    HAXX::quaternion<double> dotc = DOTC(len,&X[0],stride,&Y[0],stride);

    HAXX::quaternion<double> dotus(0.,0.,0.,0.), dotcs(0.,0.,0.,0.);
    for(auto j = 0; j < len*stride; j+= stride) {
      dotus = dotus + X[j] * Y[j];
      dotcs = dotcs + HAXX::conj(X[j]) * Y[j];
    }

    BOOST_CHECK_EQUAL(dotu,dotus);
    BOOST_CHECK_EQUAL(dotc,dotcs);
  }

  // FIXME: Need a test for when strides are not equal
}
