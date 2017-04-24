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

#ifndef _HAXX_UT_BUTF_NINCLUDED
  #include <boost/test/included/unit_test.hpp>
#else
  #include <boost/test/unit_test.hpp>
#endif

#include <boost/iterator/counting_iterator.hpp>

#include "haxx.hpp"
#include "hblas/hblas1_impl.hpp"
#include "hblas/hblas2_impl.hpp"

#include <random>
#include <iterator>
#include <iostream>
#include <limits>

// Length constants
#define HBLAS1_VECLEN 500
#define HBLAS2_MATLEN ( (HBLAS1_VECLEN) *  (HBLAS1_VECLEN) )
#define HBLAS1_RAND_MIN -20
#define HBLAS1_RAND_MAX 54

// Setup Random Number generator
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(HBLAS1_RAND_MIN,HBLAS1_RAND_MAX);

// Index list for HBLAS1 UT conformation
std::vector<int> indx(boost::counting_iterator<int>(0),
  boost::counting_iterator<int>(HBLAS1_VECLEN));

// Strides to be tested
std::vector<size_t> strides = {1,2,3,5,9};






BOOST_AUTO_TEST_CASE(hblas2_gemv)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    X(HBLAS1_VECLEN), Y(HBLAS1_VECLEN), A(HBLAS2_MATLEN);

  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : Y)
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> YC(Y);

  HAXX::quaternion<double> 
    ALPHA(dis(gen),dis(gen),dis(gen),dis(gen)), 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas2_gemv will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&X[0],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*YC[i]) / HAXX::norm(Y[i]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas2_geru)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    X(HBLAS1_VECLEN), Y(HBLAS1_VECLEN), A(HBLAS2_MATLEN);

  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : Y)
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> AC(A);

  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas2_geru will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  
  HBLAS_GERU(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {

    BOOST_CHECK_EQUAL(
      ALPHA*X[i]*Y[j] + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
      A[RANK2_INDX(i,j,HBLAS1_VECLEN)]
    );

  }
}

BOOST_AUTO_TEST_CASE(hblas2_gerc)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    X(HBLAS1_VECLEN), Y(HBLAS1_VECLEN), A(HBLAS2_MATLEN);

  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : Y)
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> AC(A);

  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas2_geru will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  
  HBLAS_GERC(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {

    BOOST_CHECK_EQUAL(
      ALPHA*X[i]*HAXX::conj(Y[j]) + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
      A[RANK2_INDX(i,j,HBLAS1_VECLEN)]
    );

  }
}
