/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */


#include "haxx_ut.hpp"
#include "hblas/hblas1.hpp"

BOOST_AUTO_TEST_SUITE(HBLAS1)

BOOST_AUTO_TEST_SUITE(HBLAS1V)

BOOST_AUTO_TEST_CASE(hblas1_swapv)
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
    HBLAS_SWAPV(len,&X[0],stride,&Y[0],stride);
  
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

BOOST_AUTO_TEST_CASE(hblas1_scalv)
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

  auto chk_SCALV = [&](int stride, int len, char SIDE, char FIELD) {
    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
          if(indx % stride == 0 and indx != len*stride ) {
            if(FIELD == 'R' or FIELD == 'r')
              return CMP_Q(tmpX[indx],X[indx]*rAlpha);
            else if(FIELD == 'C' or FIELD == 'c') {
              if(SIDE == 'R' or SIDE == 'r')
                return CMP_Q(tmpX[indx],X[indx]*cAlpha);
              else 
                return CMP_Q(tmpX[indx],cAlpha*X[indx]);
            } else {
              if(SIDE == 'R' or SIDE == 'r')
                return CMP_Q(tmpX[indx],X[indx]*hAlpha);
              else 
                return CMP_Q(tmpX[indx],hAlpha*X[indx]);
            }
          } else return (tmpX[indx] == X[indx]);
        }
      )
    );
  };

  for(auto stride : strides) {

    auto len = HBLAS1_VECLEN/stride; 

    // Right scaling by real alpha
    std::copy(X.begin(),X.end(),tmpX.begin());
    HBLAS_SCALV('R',len,rAlpha,&tmpX[0],stride);
    chk_SCALV(stride,len,'R','R');

    // Left scaling by real alpha
    std::copy(X.begin(),X.end(),tmpX.begin());
    HBLAS_SCALV('L',len,rAlpha,&tmpX[0],stride);
    chk_SCALV(stride,len,'L','R');

    // Right scaling by complex alpha
    std::copy(X.begin(),X.end(),tmpX.begin());
    HBLAS_SCALV('R',len,cAlpha,&tmpX[0],stride);
    chk_SCALV(stride,len,'R','C');

    // Left scaling by complex alpha
    std::copy(X.begin(),X.end(),tmpX.begin());
    HBLAS_SCALV('L',len,cAlpha,&tmpX[0],stride);
    chk_SCALV(stride,len,'L','C');

    // Right scaling by quaternion alpha
    std::copy(X.begin(),X.end(),tmpX.begin());
    HBLAS_SCALV('R',len,hAlpha,&tmpX[0],stride);
    chk_SCALV(stride,len,'R','H');

    // Left scaling by quaternion alpha
    std::copy(X.begin(),X.end(),tmpX.begin());
    HBLAS_SCALV('L',len,hAlpha,&tmpX[0],stride);
    chk_SCALV(stride,len,'L','H');

  }

};

BOOST_AUTO_TEST_CASE(hblas1_copyv)
{
  // Random Quaternion vector
  std::vector<HAXX::quaternion<double>> X(HBLAS1_VECLEN), Y(HBLAS1_VECLEN);
  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  // Check swap when both strides are equal
  for(int stride : strides) { 

    std::fill(Y.begin(),Y.end(),HAXX::quaternion<double>(0));

    int len = HBLAS1_VECLEN/stride;
    HBLAS_COPYV(len,&X[0],stride,&Y[0],stride);
  
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

BOOST_AUTO_TEST_CASE(hblas1_axpyv)
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
  

  auto chk_AXPYV = [&](int stride, int len, char SIDE, char FIELD) {
    BOOST_CHECK( 
      std::all_of(indx.begin(),indx.end(),
        [&](int indx) {
          if(indx % stride == 0 and indx != len*stride ) {
            if(FIELD == 'R' or FIELD == 'r')
              return CMP_Q(Y[indx],YC[indx] + X[indx]*rAlpha);
            else if(FIELD == 'C' or FIELD == 'c') {
              if(SIDE == 'R' or SIDE == 'r')
                return CMP_Q(Y[indx],YC[indx] + X[indx]*cAlpha);
              else 
                return CMP_Q(Y[indx],YC[indx] + cAlpha*X[indx]);
            } else {
              if(SIDE == 'R' or SIDE == 'r')
                return CMP_Q(Y[indx],YC[indx] + X[indx]*hAlpha);
              else 
                return CMP_Q(Y[indx],YC[indx] + hAlpha*X[indx]);
            }
          } else return (Y[indx] == YC[indx]);
        }
      )
    );
  };

  // Both strides equal
  for(auto stride : strides) {
    auto len = HBLAS1_VECLEN/stride; 

    // Right scaling by real alpha
    std::copy(YC.begin(),YC.end(),Y.begin());
    HBLAS_AXPYV('R',len,rAlpha,&X[0],stride,&Y[0],stride);
    chk_AXPYV(stride,len,'R','R');

    // Left scaling by real alpha
    std::copy(YC.begin(),YC.end(),Y.begin());
    HBLAS_AXPYV('L',len,rAlpha,&X[0],stride,&Y[0],stride);
    chk_AXPYV(stride,len,'L','R');


    // Right scaling by complex alpha
    std::copy(YC.begin(),YC.end(),Y.begin());
    HBLAS_AXPYV('R',len,cAlpha,&X[0],stride,&Y[0],stride);
    chk_AXPYV(stride,len,'R','C');

    // Left scaling by complex alpha
    std::copy(YC.begin(),YC.end(),Y.begin());
    HBLAS_AXPYV('L',len,cAlpha,&X[0],stride,&Y[0],stride);
    chk_AXPYV(stride,len,'L','C');


    // Right scaling by quaternion alpha
    std::copy(YC.begin(),YC.end(),Y.begin());
    HBLAS_AXPYV('R',len,hAlpha,&X[0],stride,&Y[0],stride);
    chk_AXPYV(stride,len,'R','H');

    // Left scaling by quaternion alpha
    std::copy(YC.begin(),YC.end(),Y.begin());
    HBLAS_AXPYV('L',len,hAlpha,&X[0],stride,&Y[0],stride);
    chk_AXPYV(stride,len,'L','H');
  }


  // FIXME: Need a test for when strides are not equal

};

BOOST_AUTO_TEST_CASE(hblas1_dotv)
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
    HAXX::quaternion<double> dotu = HBLAS_DOTUV(len,&X[0],stride,&Y[0],stride);
    HAXX::quaternion<double> dotc = HBLAS_DOTCV(len,&X[0],stride,&Y[0],stride);

    HAXX::quaternion<double> dotus(0.,0.,0.,0.), dotcs(0.,0.,0.,0.);
    for(auto j = 0; j < len*stride; j+= stride) {
      dotus = dotus + X[j] * Y[j];
      dotcs = dotcs + HAXX::conj(X[j]) * Y[j];
    }

    BOOST_CHECK_MESSAGE(CMP_Q(dotu,dotus), 
      "\n stride = " << stride << ", HBLAS_DOTUV = " << dotu 
                  << ", REFERENCE = " << dotus);
    BOOST_CHECK_MESSAGE(CMP_Q(dotc,dotcs),
      "\n stride = " << stride << ", HBLAS_DOTCV = " << dotc 
                  << ", REFERENCE = " << dotcs);
  }

  // FIXME: Need a test for when strides are not equal
}

// HBLAS1V
BOOST_AUTO_TEST_SUITE_END()

// HBLAS1
BOOST_AUTO_TEST_SUITE_END()
