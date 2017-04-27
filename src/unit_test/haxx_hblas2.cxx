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
#include "hblas/hblas2_impl.hpp"



BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NQQ_LDSame)
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

  std::cout << "hblas2_gemv_square_NQQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TQQ_LDSame)
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

  std::cout << "hblas2_gemv_square_TQQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CQQ_LDSame)
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

  std::cout << "hblas2_gemv_square_CQQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}



BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NRQ_LDSame)
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
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  double ALPHA(dis(gen));

  std::cout << "hblas2_gemv_square_NRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TRQ_LDSame)
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
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  double ALPHA(dis(gen));

  std::cout << "hblas2_gemv_square_TRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CRQ_LDSame)
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
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  double ALPHA(dis(gen));

  std::cout << "hblas2_gemv_square_CRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NCQ_LDSame)
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
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas2_gemv_square_NCQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TCQ_LDSame)
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
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas2_gemv_square_TCQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CCQ_LDSame)
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
    BETA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas2_gemv_square_CCQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NRR_LDSame)
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

  double BETA(dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas2_gemv_square_NRR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;

  
  HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TRR_LDSame)
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

  double BETA(dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas2_gemv_square_TRR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CRR_LDSame)
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

  double BETA(dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas2_gemv_square_CRR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NCR_LDSame)
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

  double BETA(dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas2_gemv_square_NCR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;

  
  HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TCR_LDSame)
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

  double BETA(dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas2_gemv_square_TCR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CCR_LDSame)
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

  double BETA(dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas2_gemv_square_CCR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NQR_LDSame)
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

  double BETA(dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas2_gemv_square_NQR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;

  
  HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(i,0,HBLAS1_VECLEN)],
        HBLAS1_VECLEN,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TQR_LDSame)
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

  double BETA(dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas2_gemv_square_TQR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CQR_LDSame)
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

  double BETA(dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas2_gemv_square_CQR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }
}





BOOST_AUTO_TEST_CASE(hblas2_geru_square_Q)
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

  std::cout << "hblas2_geru_square_Q will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  
  HBLAS_GERU(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {

    BOOST_CHECK(
      CMP_Q(
        ALPHA*X[i]*Y[j] + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
        A[RANK2_INDX(i,j,HBLAS1_VECLEN)]
      )
    );

  }
}

BOOST_AUTO_TEST_CASE(hblas2_geru_square_R)
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

  double ALPHA(dis(gen));

  std::cout << "hblas2_geru_square_R will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  
  HBLAS_GERU(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {

    BOOST_CHECK(
      CMP_Q(
        ALPHA*X[i]*Y[j] + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
        A[RANK2_INDX(i,j,HBLAS1_VECLEN)]
      )
    );

  }
}

BOOST_AUTO_TEST_CASE(hblas2_geru_square_C)
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

  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas2_geru_square_C will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  
  HBLAS_GERU(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {

    BOOST_CHECK(
      CMP_Q(
        ALPHA*X[i]*Y[j] + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
        A[RANK2_INDX(i,j,HBLAS1_VECLEN)]
      )
    );

  }
}

BOOST_AUTO_TEST_CASE(hblas2_gerc_square_Q)
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

  std::cout << "hblas2_gerc_square_Q will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  
  HBLAS_GERC(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {

    BOOST_CHECK(
      CMP_Q(
        ALPHA*X[i]*conj(Y[j]) + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
        A[RANK2_INDX(i,j,HBLAS1_VECLEN)]
      )
    );

  }
}

BOOST_AUTO_TEST_CASE(hblas2_gerc_square_R)
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

  double ALPHA(dis(gen));

  std::cout << "hblas2_gerc_square_R will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  
  HBLAS_GERC(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {

    BOOST_CHECK(
      CMP_Q(
        ALPHA*X[i]*conj(Y[j]) + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
        A[RANK2_INDX(i,j,HBLAS1_VECLEN)]
      )
    );

  }
}

BOOST_AUTO_TEST_CASE(hblas2_gerc_square_C)
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

  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas2_gerc_square_C will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  
  HBLAS_GERC(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {

    BOOST_CHECK(
      CMP_Q(
        ALPHA*X[i]*conj(Y[j]) + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
        A[RANK2_INDX(i,j,HBLAS1_VECLEN)]
      )
    );

  }
}
