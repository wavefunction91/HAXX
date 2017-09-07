/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */


#include "haxx.hpp"
#include <random>
#include <iterator>
#include <iostream>
#include <limits>
#include <chrono>
#include "hblas/hblas3_def.hpp"
#include "hblas/hblas_util_impl.hpp"

extern "C" {
  void zgemm_(const char*, const char*, const int*, const int*, 
    const int*, const std::complex<double>*, const std::complex<double>*, const int*, 
    const std::complex<double>*, const int*, const std::complex<double>*, const std::complex<double>*,
    const int*);
};

//#define GEMM_LEN 2000
#define HBLAS1_RAND_MIN -20
#define HBLAS1_RAND_MAX 54

// Setup Random Number generator
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(HBLAS1_RAND_MIN,HBLAS1_RAND_MAX);

int main() {

  for(int GEMM_LEN = 500; GEMM_LEN <= 10000; GEMM_LEN += 500) {
  std::vector<HAXX::quaternion<double>> 
    A(GEMM_LEN*GEMM_LEN), B(GEMM_LEN*GEMM_LEN), C(GEMM_LEN*GEMM_LEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  
  std::vector<std::complex<double>> 
    AC(2*GEMM_LEN*2*GEMM_LEN), BC(2*GEMM_LEN*2*GEMM_LEN), CC(2*GEMM_LEN*2*GEMM_LEN);

  int x2x = 2*GEMM_LEN;
  int gl  = GEMM_LEN;
  hzexp_(&gl,&gl,&A[0],&gl,&AC[0],&x2x);
  hzexp_(&gl,&gl,&B[0],&gl,&BC[0],&x2x);
  hzexp_(&gl,&gl,&C[0],&gl,&CC[0],&x2x);


  char TRANSA = 'N', TRANSB = 'N';

  std::complex<double> ALPHA(dis(gen),dis(gen)), BETA(dis(gen),dis(gen));
  
  if(GEMM_LEN == 500)
  HBLAS_GEMM(TRANSA,TRANSB,GEMM_LEN,GEMM_LEN,GEMM_LEN,ALPHA,&A[0],GEMM_LEN,
    &B[0],GEMM_LEN,BETA,&C[0],GEMM_LEN);

  auto hgemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM(TRANSA,TRANSB,GEMM_LEN,GEMM_LEN,GEMM_LEN,ALPHA,&A[0],GEMM_LEN,
    &B[0],GEMM_LEN,BETA,&C[0],GEMM_LEN);
  auto hgemmEnd = std::chrono::high_resolution_clock::now();


  auto zgemmStart = std::chrono::high_resolution_clock::now();
  zgemm_(&TRANSA,&TRANSB,&x2x,&x2x,&x2x,&ALPHA,&AC[0],&x2x,
    &BC[0],&x2x,&BETA,&CC[0],&x2x);
  auto zgemmEnd = std::chrono::high_resolution_clock::now();
  
  std::chrono::duration<double> hgemmDur = hgemmEnd - hgemmStart;
  std::chrono::duration<double> zgemmDur = zgemmEnd - zgemmStart;

  std::cout << "HGEMM " << GEMM_LEN << " " << hgemmDur.count() << std::endl;
  std::cout << "ZGEMM " << GEMM_LEN << " " << zgemmDur.count() << std::endl;
  }
}
