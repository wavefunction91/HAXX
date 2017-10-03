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
#include <iomanip>
#include <limits>
#include <chrono>
#include "hblas/hblas3_def.hpp"
#include "hblas/hblas_util.hpp"

extern "C" {

  void zgemm_(const char*, const char*, const int*, const int*, 
    const int*, const std::complex<double>*, const std::complex<double>*, 
    const int*, const std::complex<double>*, const int*, 
    const std::complex<double>*, const std::complex<double>*, const int*);

};

#define HBLAS1_RAND_MIN -20
#define HBLAS1_RAND_MAX 54

// Setup Random Number generator
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(HBLAS1_RAND_MIN,HBLAS1_RAND_MAX);


//#define _DO_COMPLEX
//#define _DO_FORTRAN


void outTime(std::string name, size_t LEN, size_t FLOPS, double count) {

  std::cout << std::setw(15) << std::left  << name;
  std::cout << std::setw(15) << std::left << LEN;

  std::cout << std::setprecision(8);

  std::cout << std::setw(15) << std::right << count;
  std::cout << std::setw(15) << std::right << FLOPS/count/1.e9;

  std::cout << std::endl;



};

int main() {

  for(int GEMM_LEN = 500; GEMM_LEN <= 2500; GEMM_LEN += 500) {
  std::vector<HAXX::quaternion<double>> 
    A(GEMM_LEN*GEMM_LEN), B(GEMM_LEN*GEMM_LEN), C(GEMM_LEN*GEMM_LEN);

  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : B) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : C) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  
#ifdef _DO_COMPLEX
  std::vector<std::complex<double>> 
    AC(2*GEMM_LEN*2*GEMM_LEN), BC(2*GEMM_LEN*2*GEMM_LEN), CC(2*GEMM_LEN*2*GEMM_LEN);

  int x2x = 2*GEMM_LEN;
  int gl  = GEMM_LEN;
  HBLAS_COMPLEX_EXPAND('S',GEMM_LEN,GEMM_LEN,&A[0],GEMM_LEN,&AC[0],2*GEMM_LEN);
  HBLAS_COMPLEX_EXPAND('S',GEMM_LEN,GEMM_LEN,&B[0],GEMM_LEN,&BC[0],2*GEMM_LEN);
  HBLAS_COMPLEX_EXPAND('S',GEMM_LEN,GEMM_LEN,&C[0],GEMM_LEN,&CC[0],2*GEMM_LEN);
#endif


  char TRANSA = 'N', TRANSB = 'N';

  std::complex<double> ALPHA(dis(gen),dis(gen)), BETA(dis(gen),dis(gen));
  
  if(GEMM_LEN == 500)
  HBLAS_GEMM(TRANSA,TRANSB,GEMM_LEN,GEMM_LEN,GEMM_LEN,ALPHA,&A[0],GEMM_LEN,
    &B[0],GEMM_LEN,BETA,&C[0],GEMM_LEN);

  auto hgemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM(TRANSA,TRANSB,GEMM_LEN,GEMM_LEN,GEMM_LEN,ALPHA,&A[0],GEMM_LEN,
    &B[0],GEMM_LEN,BETA,&C[0],GEMM_LEN);
  auto hgemmEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> hgemmDur = hgemmEnd - hgemmStart;

  outTime("HGEMM",GEMM_LEN,32.*GEMM_LEN*GEMM_LEN*GEMM_LEN,hgemmDur.count());

#ifdef _DO_FORTRAN
  auto fortranStart = std::chrono::high_resolution_clock::now();
  hgemmzz_(&TRANSA,&TRANSB,&gl,&gl,&gl,&ALPHA,&A[0],&gl,&B[0],&gl,&BETA,&C[0],&gl);
  auto fortranEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> fortranDur = fortranEnd - fortranStart;

  outTime("HGEMM_FORTRAN",GEMM_LEN,32.*GEMM_LEN*GEMM_LEN*GEMM_LEN,
    fortranDur.count());
#endif

#ifdef _DO_COMPLEX
  auto zgemmStart = std::chrono::high_resolution_clock::now();
  zgemm_(&TRANSA,&TRANSB,&x2x,&x2x,&x2x,&ALPHA,&AC[0],&x2x,
    &BC[0],&x2x,&BETA,&CC[0],&x2x);
  auto zgemmEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> zgemmDur = zgemmEnd - zgemmStart;

  outTime("ZGEMM",GEMM_LEN,64.*GEMM_LEN*GEMM_LEN*GEMM_LEN,zgemmDur.count());
#endif
  

  }
}
