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

#include "haxx_ut.hpp"
#include "hblas/hblas1_impl.hpp"
#include "hblas/hblas2_impl.hpp"
#include "hblas/hblas3_impl.hpp"



BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNRR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_NNRR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NNRR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNRR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_TNRR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TNRR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNRR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_CNRR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CNRR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTRR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_NTRR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NTRR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCRR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_NCRR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NCRR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 
  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTRR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_TTRR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TTRR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCRR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_TCRR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TCRR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTRR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_CTRR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CTRR_LDSame GEMM took " << gemmDur.count() << " s\n";


  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCRR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_CCRR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CCRR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNCR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NNCR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NNCR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNCR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_TNCR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TNCR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNCR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_CNCR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CNCR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTCR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NTCR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NTCR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCCR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NCCR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NCCR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 
  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTCR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_TTCR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TTCR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCCR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_TCCR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TCCR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTCR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_CTCR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CTCR_LDSame GEMM took " << gemmDur.count() << " s\n";


  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCCR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_CCCR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CCCR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNQR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NNQR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NNQR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNQR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_TNQR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TNQR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNQR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_CNQR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CNQR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTQR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NTQR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NTQR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCQR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NCQR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NCQR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 
  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTQR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_TTQR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TTQR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCQR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_TCQR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TCQR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTQR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_CTQR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CTQR_LDSame GEMM took " << gemmDur.count() << " s\n";


  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCQR_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  double BETA(dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_CCQR_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CCQR_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}



BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNRC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_NNRC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NNRC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNRC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_TNRC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TNRC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNRC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_CNRC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CNRC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTRC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_NTRC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NTRC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCRC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_NCRC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NCRC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 
  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTRC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_TTRC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TTRC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCRC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_TCRC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TCRC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTRC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_CTRC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CTRC_LDSame GEMM took " << gemmDur.count() << " s\n";


  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCRC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_CCRC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CCRC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNCC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NNCC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NNCC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNCC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_TNCC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TNCC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNCC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_CNCC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CNCC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTCC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NTCC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NTCC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCCC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NCCC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NCCC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 
  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTCC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_TTCC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TTCC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCCC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_TCCC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TCCC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTCC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_CTCC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CTCC_LDSame GEMM took " << gemmDur.count() << " s\n";


  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCCC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_CCCC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CCCC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNQC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NNQC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NNQC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNQC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_TNQC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TNQC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNQC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_CNQC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CNQC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTQC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NTQC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NTQC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCQC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NCQC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NCQC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 
  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTQC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_TTQC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TTQC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCQC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_TCQC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TCQC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTQC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_CTQC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CTQC_LDSame GEMM took " << gemmDur.count() << " s\n";


  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCQC_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  std::complex<double> BETA(dis(gen),dis(gen));
  HAXX::quaternion<double> ALPHA(dis(gen),dis(gen),dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_CCQC_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CCQC_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_NNRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NNRQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_TNRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TNRQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_CNRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CNRQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_NTRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NTRQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_NCRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NCRQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 
  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_TTRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TTRQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_TCRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TCRQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_CTRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CTRQ_LDSame GEMM took " << gemmDur.count() << " s\n";


  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  double ALPHA(dis(gen));

  std::cout << "hblas3_gemm_square_CCRQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CCRQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNZQ_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NNZQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NNZQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNZQ_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_TNZQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TNZQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNZQ_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_CNZQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CNZQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j*HBLAS1_VECLEN],1,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTZQ_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NTZQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NTZQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCZQ_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_NCZQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('N','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_NCZQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 
  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('N',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTZQ_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_TTZQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TTZQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }

}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCZQ_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_TCZQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('T','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_TCZQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('T',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTZQ_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_CTZQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','T',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CTZQ_LDSame GEMM took " << gemmDur.count() << " s\n";


  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCZQ_LDSame)
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
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  HAXX::quaternion<double> BETA(dis(gen),dis(gen),dis(gen),dis(gen));
  std::complex<double> ALPHA(dis(gen),dis(gen));

  std::cout << "hblas3_gemm_square_CCZQ_LDSame will use " << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;
  
  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM('C','C',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << "hblas3_gemm_square_CCZQ_LDSame GEMM took " << gemmDur.count() << " s\n";

  for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HBLAS_GEMV('C',HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[j],HBLAS1_VECLEN,0.,&SCR[0],1);
    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
}
