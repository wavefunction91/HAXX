/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef _HAXX_UT_BUTF_NINCLUDED
  #include <boost/test/included/unit_test.hpp>
#else
  #include <boost/test/unit_test.hpp>
#endif

#include <boost/iterator/counting_iterator.hpp>

#include "haxx.hpp"

#include <random>
#include <iterator>
#include <iostream>
#include <limits>
#include <chrono>

// Length constants
#define HBLAS1_VECLEN 500
#define HBLAS2_MATLEN (HBLAS1_VECLEN) * (HBLAS1_VECLEN)
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

