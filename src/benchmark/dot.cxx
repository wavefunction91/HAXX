/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#include "haxx.hpp"
#include "hblas/hblas1.hpp"

#include <random>
#include <iterator>
#include <iostream>
#include <limits>
#include <chrono>


#define HBLAS1_RAND_MIN -20
#define HBLAS1_RAND_MAX 54

// Setup Random Number generator
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(HBLAS1_RAND_MIN,HBLAS1_RAND_MAX);

#define DOT_LEN_MAX 10000000
#define DOT_LEN_START 100
#define NTEST 10
#define NREP 10

#define USE_ALIGNED_MEM

using namespace HAXX;

int main() {


#ifdef USE_ALIGNED_MEM

  size_t alignment = 32;
  size_t nAlloc = (size_t(DOT_LEN_MAX) * sizeof(quaternion<double>));
  nAlloc = nAlloc + (nAlloc % 32);

  quaternion<double> *X = 
    (quaternion<double>*)aligned_alloc(alignment, nAlloc);
  quaternion<double> *Y = 
    (quaternion<double>*)aligned_alloc(alignment, nAlloc);

  for(auto i = 0ul; i < DOT_LEN_MAX; i++) {
    X[i] = quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
    Y[i] = quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  }

#else

  std::vector<quaternion<double>> X(DOT_LEN_MAX), Y(DOT_LEN_MAX);
  for(auto &x : X) x = quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : Y) x = quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

#endif


  for(auto N = DOT_LEN_START; N <= DOT_LEN_MAX; 
    N += (DOT_LEN_MAX - DOT_LEN_START)/NTEST ) {

    // Fix cache
    if( N == DOT_LEN_START ) HBLAS_DOTUV(N,&X[0],1,&Y[0],1);

    std::chrono::duration<double> dur_opt(0.);  
    std::chrono::duration<double> dur_for(0.);  

    quaternion<double> res;


    int INCX = 1;
    for(auto rep = 0; rep < NREP; rep++) {
      auto start = std::chrono::high_resolution_clock::now();
      hdotu_(&res,&N,&X[0],&INCX,&Y[0],&INCX);
      auto end = std::chrono::high_resolution_clock::now();
      dur_for += end - start;  
    }

    for(auto rep = 0; rep < NREP; rep++) {
      auto start = std::chrono::high_resolution_clock::now();
      res  = HBLAS_DOTUV(N,&X[0],1,&Y[0],1);
      auto end = std::chrono::high_resolution_clock::now();
      dur_opt += end - start;  
    }


    dur_opt /= NREP;
    dur_for /= NREP;

    std::cout << "N = " << N << ", SIMD = " << dur_opt.count() 
              << ", FORTRAN = " << dur_for.count() 
              << ", % = " <<  (dur_for.count() - dur_opt.count()) / dur_for.count() 
              << std::endl;

  }

#ifdef USE_ALIGNED_MEM
  free(X); free(Y);
#endif

};
