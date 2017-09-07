/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#include "haxx_ut.hpp"
#include "hblas/hblas_util.hpp"
#include "hblas/hblas3_def.hpp"

BOOST_AUTO_TEST_SUITE(UTIL)

void ComplexExpandTest(char ORDER) {

  std::vector<HAXX::quaternion<double>> A(HBLAS2_MATLEN);
  std::vector<HAXX::quaternion<double>> B(HBLAS2_MATLEN);
  std::vector<HAXX::quaternion<double>> C(HBLAS2_MATLEN);
  
  std::vector<std::complex<double>> AC(4*HBLAS2_MATLEN);
  std::vector<std::complex<double>> BC(4*HBLAS2_MATLEN);
  std::vector<std::complex<double>> CC(4*HBLAS2_MATLEN);
  std::vector<std::complex<double>> PC(4*HBLAS2_MATLEN);

  for(auto &x : A) x = genRandom<HAXX::quaternion<double>>();
  for(auto &x : B) x = genRandom<HAXX::quaternion<double>>();

  HBLAS_COMPLEX_EXPAND(ORDER,HBLAS1_VECLEN,HBLAS1_VECLEN,&A[0],HBLAS1_VECLEN,
    &AC[0],2*HBLAS1_VECLEN);
  HBLAS_COMPLEX_EXPAND(ORDER,HBLAS1_VECLEN,HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,
    &BC[0],2*HBLAS1_VECLEN);


  // Quaternion multiplication
  HBLAS_GEMM('N','N', HBLAS1_VECLEN, HBLAS1_VECLEN, HBLAS1_VECLEN,
    1., &A[0], HBLAS1_VECLEN, &B[0], HBLAS1_VECLEN, 0., &C[0],
    HBLAS1_VECLEN);
  
  // Complex multiplication (FIXME: Should use ZGEMM here)
  for(auto i = 0; i < 2*HBLAS1_VECLEN; i++)
  for(auto j = 0; j < 2*HBLAS1_VECLEN; j++) {
    CC[i + j*2*HBLAS1_VECLEN]  = 0.;
  for(auto k = 0; k < 2*HBLAS1_VECLEN; k++)
    CC[i + j*2*HBLAS1_VECLEN]  += 
      AC[i + k*2*HBLAS1_VECLEN] * BC[k + j*2*HBLAS1_VECLEN];
  }

  // Expand product to complex
  HBLAS_COMPLEX_EXPAND(ORDER,HBLAS1_VECLEN,HBLAS1_VECLEN,&C[0],HBLAS1_VECLEN,
    &PC[0],2*HBLAS1_VECLEN);
  

  for(auto i = 0; i < 4*HBLAS2_MATLEN; i++)
    BOOST_CHECK( std::norm(CC[i] / PC[i] - 1.) < COMPARE_TOL );

};

void ComplexContractTest(char ORDER, char UPLO) {

  std::vector<HAXX::quaternion<double>> A(HBLAS2_MATLEN);
  std::vector<HAXX::quaternion<double>> B(HBLAS2_MATLEN);
  
  std::vector<std::complex<double>> AC(4*HBLAS2_MATLEN);

  for(auto &x : A) x = genRandom<HAXX::quaternion<double>>();

  // A -> AC
  HBLAS_COMPLEX_EXPAND(ORDER,HBLAS1_VECLEN,HBLAS1_VECLEN,&A[0],HBLAS1_VECLEN,
    &AC[0],2*HBLAS1_VECLEN);

  // AC -> B
  HBLAS_COMPLEX_CONTRACT(ORDER,UPLO,HBLAS1_VECLEN,HBLAS1_VECLEN,&B[0],
    HBLAS1_VECLEN,&AC[0],2*HBLAS1_VECLEN);


  // Compare A and B
  for(auto i = 0; i < HBLAS2_MATLEN; i++)
    BOOST_CHECK( CMP_Q(A[i],B[i]) );
}


BOOST_AUTO_TEST_CASE(Complex_Expand1) { ComplexExpandTest('F'); };
BOOST_AUTO_TEST_CASE(Complex_Expand2) { ComplexExpandTest('S'); };
BOOST_AUTO_TEST_CASE(Complex_Contract1) { ComplexContractTest('F','U'); };
BOOST_AUTO_TEST_CASE(Complex_Contract2) { ComplexContractTest('S','U'); };
BOOST_AUTO_TEST_CASE(Complex_Contract3) { ComplexContractTest('F','L'); };
BOOST_AUTO_TEST_CASE(Complex_Contract4) { ComplexContractTest('S','L'); };

BOOST_AUTO_TEST_SUITE_END()
