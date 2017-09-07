/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */


#ifdef BOOST_TEST_MODULE
  #undef BOOST_TEST_MODULE
#endif

#define BOOST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include "haxx.hpp"

#include <random>
#include <iterator>
#include <iostream>
#include <limits>
#include <chrono>

// Length constants
#define HBLAS1_VECLEN 100
#define HBLAS2_MATLEN (HBLAS1_VECLEN) * (HBLAS1_VECLEN)
#define HBLAS1_RAND_MIN -20
#define HBLAS1_RAND_MAX 54

// Setup Random Number generator
static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<> dis(HBLAS1_RAND_MIN,HBLAS1_RAND_MAX);

template <typename _F> _F genRandom();
template<> inline double genRandom<double>(){ return double(dis(gen)); }
template<> inline std::complex<double> genRandom<std::complex<double>>(){ 
  return std::complex<double>(dis(gen),dis(gen)); 
}
template<> inline HAXX::quaternion<double> genRandom<HAXX::quaternion<double>>(){ 
  return HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen)); 
}



// Index list for HBLAS1 UT conformation
static std::vector<int> indx(boost::counting_iterator<int>(0),
  boost::counting_iterator<int>(HBLAS1_VECLEN));

// Strides to be tested
static std::vector<size_t> strides = {1,2,3,5,9};

#define COMPARE_TOL 1e-12
//#define CMP_Q(a,b) (HAXX::norm(a)/HAXX::norm(b) - 1. < COMPARE_TOL)
#define CMP_Q(a,b) ( HAXX::norm(((a) * HAXX::inv(b))- 1.) < COMPARE_TOL )
