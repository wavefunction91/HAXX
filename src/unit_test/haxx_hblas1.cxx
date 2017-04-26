/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HAXX_BLAS1

#include "haxx_ut.hpp"
#include "hblas/hblas1_impl.hpp"


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

  // Check swap when both strides are equal
  for(int stride : strides) { 

    std::copy(XC.begin(),XC.end(),X.begin());
    std::copy(YC.begin(),YC.end(),Y.begin());

    int len = HBLAS1_VECLEN/stride;
    HBLAS_SWAP(len,&X[0],stride,&Y[0],stride);
  
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

  for(auto stride : strides) {


    // Right scaling by real alpha
    std::copy(X.begin(),X.end(),tmpX.begin());

    auto len = HBLAS1_VECLEN/stride; 
    HBLAS_SCAL('R',len,rAlpha,&tmpX[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (HAXX::norm(tmpX[indx])/HAXX::norm(X[indx]*rAlpha)  - 1. < 1e-12);
	  else return (tmpX[indx] == X[indx]);
        }
      )
    );


    // Left scaling by real alpha
    std::copy(X.begin(),X.end(),tmpX.begin());

    HBLAS_SCAL('L',len,rAlpha,&tmpX[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (HAXX::norm(tmpX[indx])/HAXX::norm(X[indx]*rAlpha)  - 1. < 1e-12);
	  else return (tmpX[indx] == X[indx]);
        }
      )
    );

    // Right scaling by complex alpha
    std::copy(X.begin(),X.end(),tmpX.begin());

    HBLAS_SCAL('R',len,cAlpha,&tmpX[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (HAXX::norm(tmpX[indx])/HAXX::norm(X[indx]*cAlpha)  - 1. < 1e-12);
	  else return (tmpX[indx] == X[indx]);
        }
      )
    );


    // Left scaling by complex alpha
    std::copy(X.begin(),X.end(),tmpX.begin());

    HBLAS_SCAL('L',len,cAlpha,&tmpX[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (HAXX::norm(tmpX[indx])/HAXX::norm(cAlpha*X[indx])  - 1. < 1e-12);
	  else return (tmpX[indx] == X[indx]);
        }
      )
    );

    // Right scaling by quaternion alpha
    std::copy(X.begin(),X.end(),tmpX.begin());

    HBLAS_SCAL('R',len,hAlpha,&tmpX[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (HAXX::norm(tmpX[indx])/HAXX::norm(X[indx]*hAlpha)  - 1. < 1e-12);
	  else return (tmpX[indx] == X[indx]);
        }
      )
    );


    // Left scaling by quaternion alpha
    std::copy(X.begin(),X.end(),tmpX.begin());

    HBLAS_SCAL('L',len,hAlpha,&tmpX[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (HAXX::norm(tmpX[indx])/HAXX::norm(hAlpha*X[indx])  - 1. < 1e-12);
	  else return (tmpX[indx] == X[indx]);
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

  // Check swap when both strides are equal
  for(int stride : strides) { 

    std::fill(Y.begin(),Y.end(),HAXX::quaternion<double>(0));

    int len = HBLAS1_VECLEN/stride;
    HBLAS_COPY(len,&X[0],stride,&Y[0],stride);
  
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

  std::vector<HAXX::quaternion<double>> YC(Y);
  
  // Both strides equal
  for(auto stride : strides) {
    auto len = HBLAS1_VECLEN/stride; 

    // Right scaling by real alpha
    std::copy(YC.begin(),YC.end(),Y.begin());

    HBLAS_AXPY('R',len,rAlpha,&X[0],stride,&Y[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (HAXX::norm(Y[indx])/HAXX::norm(YC[indx] + X[indx]*rAlpha) - 1. < 1e-12);
	  else return (Y[indx] == YC[indx]);
        }
      )
    );

    // Left scaling by real alpha
    std::copy(YC.begin(),YC.end(),Y.begin());

    HBLAS_AXPY('L',len,rAlpha,&X[0],stride,&Y[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (HAXX::norm(Y[indx])/HAXX::norm(YC[indx] + X[indx]*rAlpha) - 1. < 1e-12);
	  else return (Y[indx] == YC[indx]);
        }
      )
    );

    // Right scaling by complex alpha
    std::copy(YC.begin(),YC.end(),Y.begin());

    HBLAS_AXPY('R',len,cAlpha,&X[0],stride,&Y[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (HAXX::norm(Y[indx])/HAXX::norm(YC[indx] + X[indx]*cAlpha) - 1. < 1e-12);
	  else return (Y[indx] == YC[indx]);
        }
      )
    );

    // Left scaling by complex alpha
    std::copy(YC.begin(),YC.end(),Y.begin());

    HBLAS_AXPY('L',len,cAlpha,&X[0],stride,&Y[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (HAXX::norm(Y[indx])/HAXX::norm(YC[indx] + cAlpha*X[indx]) - 1. < 1e-12);
	  else return (Y[indx] == YC[indx]);
        }
      )
    );


    // Right scaling by quaternion alpha
    std::copy(YC.begin(),YC.end(),Y.begin());

    HBLAS_AXPY('R',len,hAlpha,&X[0],stride,&Y[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (HAXX::norm(Y[indx])/HAXX::norm(YC[indx] + X[indx]*hAlpha) - 1. < 1e-12);
	  else return (Y[indx] == YC[indx]);
        }
      )
    );

    // Left scaling by quaternion alpha
    std::copy(YC.begin(),YC.end(),Y.begin());

    HBLAS_AXPY('L',len,hAlpha,&X[0],stride,&Y[0],stride);

    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
	  if(indx % stride == 0 and indx != len*stride ) 
            return (HAXX::norm(Y[indx])/HAXX::norm(YC[indx] + hAlpha*X[indx]) - 1. < 1e-12);
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


  // Both strides equal
  for(int stride : strides) {
    auto len = HBLAS1_VECLEN/stride; 
    HAXX::quaternion<double> dotu = HBLAS_DOTU(len,&X[0],stride,&Y[0],stride);
    HAXX::quaternion<double> dotc = HBLAS_DOTC(len,&X[0],stride,&Y[0],stride);

    HAXX::quaternion<double> dotus(0.,0.,0.,0.), dotcs(0.,0.,0.,0.);
    for(auto j = 0; j < len*stride; j+= stride) {
      dotus = dotus + X[j] * Y[j];
      dotcs = dotcs + HAXX::conj(X[j]) * Y[j];
    }

    BOOST_CHECK_CLOSE(1.,HAXX::norm(dotu)/HAXX::norm(dotus),1e-12);
    BOOST_CHECK_CLOSE(1.,HAXX::norm(dotc)/HAXX::norm(dotcs),1e-12);
  }

  // FIXME: Need a test for when strides are not equal
}
