/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HAXX_BLAS3

#ifndef _HAXX_UT_BUTF_NINCLUDED
  #include <boost/test/included/unit_test.hpp>
#else
  #include <boost/test/unit_test.hpp>
#endif

#include <boost/iterator/counting_iterator.hpp>

#include "haxx.hpp"
#include "hblas/hblas1_impl.hpp"
#include "hblas/hblas2_impl.hpp"
#include "hblas/hblas3_impl.hpp"

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






BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNQQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    ALPHA(dis(gen),dis(gen),dis(gen),dis(gen)), 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NNQQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);


  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&B[RANK2_INDX(0,j,HBLAS1_VECLEN)],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNQQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    ALPHA(dis(gen),dis(gen),dis(gen),dis(gen)), 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_TNQQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('T','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);


  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&B[RANK2_INDX(0,j,HBLAS1_VECLEN)],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNQQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    ALPHA(dis(gen),dis(gen),dis(gen),dis(gen)), 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_CNQQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('C','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);


  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&B[RANK2_INDX(0,j,HBLAS1_VECLEN)],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTQQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    ALPHA(dis(gen),dis(gen),dis(gen),dis(gen)), 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NTQQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('N','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);


  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&B[RANK2_INDX(j,0,HBLAS1_VECLEN)],HBLAS1_VECLEN);


    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCQQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    ALPHA(dis(gen),dis(gen),dis(gen),dis(gen)), 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NCQQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('N','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&B[RANK2_INDX(j,0,HBLAS1_VECLEN)],HBLAS1_VECLEN);


    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTQQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    ALPHA(dis(gen),dis(gen),dis(gen),dis(gen)), 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_TTQQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('T','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);


  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&B[RANK2_INDX(j,0,HBLAS1_VECLEN)],HBLAS1_VECLEN);


    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCQQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    ALPHA(dis(gen),dis(gen),dis(gen),dis(gen)), 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_TCQQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('T','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&B[RANK2_INDX(j,0,HBLAS1_VECLEN)],HBLAS1_VECLEN);


    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTQQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    ALPHA(dis(gen),dis(gen),dis(gen),dis(gen)), 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_CTQQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('C','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);


  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&B[RANK2_INDX(j,0,HBLAS1_VECLEN)],HBLAS1_VECLEN);


    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCQQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    ALPHA(dis(gen),dis(gen),dis(gen),dis(gen)), 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_CCQQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('C','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&B[RANK2_INDX(j,0,HBLAS1_VECLEN)],HBLAS1_VECLEN);


    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}






BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNRQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_NNRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);


  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&B[RANK2_INDX(0,j,HBLAS1_VECLEN)],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNRQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_TNRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('T','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);


  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&B[RANK2_INDX(0,j,HBLAS1_VECLEN)],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNRQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_CNRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('C','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);


  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&B[RANK2_INDX(0,j,HBLAS1_VECLEN)],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTRQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_NTRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('N','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);


  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&B[RANK2_INDX(j,0,HBLAS1_VECLEN)],HBLAS1_VECLEN);


    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCRQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_NCRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('N','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&B[RANK2_INDX(j,0,HBLAS1_VECLEN)],HBLAS1_VECLEN);


    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTRQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_TTRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('T','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);


  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&B[RANK2_INDX(j,0,HBLAS1_VECLEN)],HBLAS1_VECLEN);


    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCRQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_TCRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('T','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&B[RANK2_INDX(j,0,HBLAS1_VECLEN)],HBLAS1_VECLEN);


    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTRQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_CTRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('C','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);


  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&B[RANK2_INDX(j,0,HBLAS1_VECLEN)],HBLAS1_VECLEN);


    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCRQ_LDSame)
{
  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), B(HBLAS2_MATLEN), C(HBLAS2_MATLEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> CC(C);

  HAXX::quaternion<double> 
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_CCRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMM('C','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&B[RANK2_INDX(j,0,HBLAS1_VECLEN)],HBLAS1_VECLEN);


    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]) / 
      HAXX::norm(C[RANK2_INDX(i,j,HBLAS1_VECLEN)]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}
