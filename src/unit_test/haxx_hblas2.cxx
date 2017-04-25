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



BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NQ_LDSame)
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

  std::cout << "hblas2_gemv_square_NQ_LDSame will use " << std::endl;
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

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TQ_LDSame)
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

  std::cout << "hblas2_gemv_square_TQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*YC[i]) / HAXX::norm(Y[i]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CQ_LDSame)
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

  std::cout << "hblas2_gemv_square_CQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*YC[i]) / HAXX::norm(Y[i]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}



BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NR_LDSame)
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

  std::cout << "hblas2_gemv_square_NR_LDSame will use " << std::endl;
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

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TR_LDSame)
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

  std::cout << "hblas2_gemv_square_TR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*YC[i]) / HAXX::norm(Y[i]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CR_LDSame)
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

  std::cout << "hblas2_gemv_square_CR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*YC[i]) / HAXX::norm(Y[i]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NC_LDSame)
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

  std::cout << "hblas2_gemv_square_NR_LDSame will use " << std::endl;
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

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TC_LDSame)
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

  std::cout << "hblas2_gemv_square_TR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*YC[i]) / HAXX::norm(Y[i]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
  }
}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CC_LDSame)
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

  std::cout << "hblas2_gemv_square_CR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX::quaternion<double>
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(0,i,HBLAS1_VECLEN)],
        1,&X[0],1);

    // FIXME: epsilon check is too tight, what is a proper criteria here
    //   in relation to machine epsilon?
    BOOST_CHECK_CLOSE(double(1.), 
      HAXX::norm(ALPHA*tmp + BETA*YC[i]) / HAXX::norm(Y[i]),
      //std::numeric_limits<double>::epsilon() * 4);
      1e-12);
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

    BOOST_CHECK_EQUAL(
      ALPHA*X[i]*Y[j] + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
      A[RANK2_INDX(i,j,HBLAS1_VECLEN)]
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

    BOOST_CHECK_CLOSE(double(1.),
      HAXX::norm(ALPHA*X[i]*Y[j] + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)])/
      HAXX::norm(A[RANK2_INDX(i,j,HBLAS1_VECLEN)]),1e-12
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

  std::cout << "hblas2_geru_square_R will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  
  HBLAS_GERU(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {

    BOOST_CHECK_CLOSE(double(1.),
      HAXX::norm(ALPHA*X[i]*Y[j] + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)])/
      HAXX::norm(A[RANK2_INDX(i,j,HBLAS1_VECLEN)]),1e-12
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

    BOOST_CHECK_EQUAL(
      ALPHA*X[i]*HAXX::conj(Y[j]) + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
      A[RANK2_INDX(i,j,HBLAS1_VECLEN)]
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

  std::cout << "hblas2_geru_square_R will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  
  HBLAS_GERC(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {

    BOOST_CHECK_CLOSE(double(1.),
      HAXX::norm(ALPHA*X[i]*HAXX::conj(Y[j]) + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)])/
      HAXX::norm(A[RANK2_INDX(i,j,HBLAS1_VECLEN)]),1e-12
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

  std::cout << "hblas2_geru_square_C will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  
  HBLAS_GERC(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {

    BOOST_CHECK_CLOSE(double(1.),
      HAXX::norm(ALPHA*X[i]*HAXX::conj(Y[j]) + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)])/
      HAXX::norm(A[RANK2_INDX(i,j,HBLAS1_VECLEN)]),1e-12
    );

  }
}
