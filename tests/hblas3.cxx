/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#include "haxx_ut.hpp"
#include "hblas/hblas1.hpp"
#include "hblas/hblas2.hpp"
#include "hblas/hblas3.hpp"

BOOST_AUTO_TEST_SUITE(HBLAS3)


template <typename _AF, typename _BF, typename _AlphaF, typename _BetaF, char _TA, char _TB>
void hblas3_gemm_square_LDSame() {

  // Random Quaternion vectors and matricies
  std::vector<_AF> A(HBLAS2_MATLEN);
  std::vector<_BF> B(HBLAS2_MATLEN);
  std::vector<HAXX::quaternion<double>> C(HBLAS2_MATLEN);

  for(auto &x : A) x = genRandom<_AF>();
  for(auto &x : B) x = genRandom<_BF>();
  for(auto &x : C) x = genRandom<HAXX::quaternion<double>>();

  std::vector<HAXX::quaternion<double>> CC(C);
  std::vector<HAXX::quaternion<double>> SCR(HBLAS1_VECLEN);

  _AlphaF ALPHA = genRandom<_AlphaF>();
  _BetaF  BETA = genRandom<_BetaF>();

  std::stringstream ss;

  ss << "hblas3_gemm_square_";
  if(std::is_same<double,_AF>::value) ss << "R";
  else if(std::is_same<std::complex<double>,_AF>::value) ss << "C";
  else ss << "Q";
  if(std::is_same<double,_BF>::value) ss << "R";
  else if(std::is_same<std::complex<double>,_BF>::value) ss << "C";
  else ss << "Q";


  ss << _TA << _TB;



  if(std::is_same<double,_AlphaF>::value) ss << "R";
  else if(std::is_same<std::complex<double>,_AlphaF>::value) ss << "C";
  else ss << "Q";
  if(std::is_same<double,_BetaF>::value) ss << "R";
  else if(std::is_same<std::complex<double>,_BetaF>::value) ss << "C";
  else ss << "Q";
  ss << "_LDSame";

  std::cout << ss.str() << " will use" << std::endl;
  std::cout << "  ALPHA = " << ALPHA << std::endl;
  std::cout << "  BETA = " << BETA << std::endl;



  auto gemmStart = std::chrono::high_resolution_clock::now();
  HBLAS_GEMM(_TA,_TB,HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
    HBLAS1_VECLEN,&B[0],HBLAS1_VECLEN,BETA,&C[0],HBLAS1_VECLEN);
  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << std::endl << ss.str() << " GEMM took " << gemmDur.count() 
    << " s\n";

  HAXX_INT BStride = 1;
  if( _TB == 'C' or _TB == 'T' ) BStride = HBLAS1_VECLEN;
  if( _TB == 'C' ) 
    for(auto &x : B) x = HAXX::conj(x); 

  for(int j = 0; j < HBLAS1_VECLEN; j++) {
    HAXX_INT BStart = j;
    if( _TB == 'N' ) BStart *= HBLAS1_VECLEN;

    HBLAS_GEMV(_TA,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],
      HBLAS1_VECLEN,&B[BStart],BStride,0.,&SCR[0],1);

    for(auto i = 0; i < HBLAS1_VECLEN; i++) {
      BOOST_CHECK_MESSAGE( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        ),
        "\n i = " << i << " j = " << j << ", HBLAS_GEMM = " 
                  << C[RANK2_INDX(i,j,HBLAS1_VECLEN)] << ", REFERENCE = "
                  << SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)]
      );
    }
  }

};



// Quaternion-Quaternion matrix-matrix multiplication

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNNRR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,double,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNTRR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,double,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNCRR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,double,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTNRR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,double,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTTRR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,double,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTCRR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,double,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCNRR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,double,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCTRR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,double,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCCRR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,double,'C','C'>();
}


BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNNCR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,double,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNTCR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,double,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNCCR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,double,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTNCR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,double,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTTCR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,double,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTCCR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,double,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCNCR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,double,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCTCR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,double,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCCCR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,double,'C','C'>();
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNNQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,double,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNTQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,double,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNCQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,double,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTNQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,double,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTTQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,double,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTCQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,double,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCNQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,double,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCTQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,double,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCCQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,double,'C','C'>();
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNNRC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,std::complex<double>,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNTRC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,std::complex<double>,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNCRC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,std::complex<double>,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTNRC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,std::complex<double>,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTTRC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,std::complex<double>,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTCRC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,std::complex<double>,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCNRC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,std::complex<double>,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCTRC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,std::complex<double>,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCCRC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,std::complex<double>,'C','C'>();
}


BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNNCC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,std::complex<double>,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNTCC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,std::complex<double>,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNCCC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,std::complex<double>,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTNCC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,std::complex<double>,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTTCC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,std::complex<double>,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTCCC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,std::complex<double>,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCNCC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,std::complex<double>,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCTCC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,std::complex<double>,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCCCC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,std::complex<double>,'C','C'>();
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNNQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNTQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNCQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTNQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTTQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTCQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCNQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCTQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCCQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,'C','C'>();
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNNRQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,HAXX::quaternion<double>,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNTRQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,HAXX::quaternion<double>,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNCRQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,HAXX::quaternion<double>,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTNRQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,HAXX::quaternion<double>,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTTRQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,HAXX::quaternion<double>,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTCRQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,HAXX::quaternion<double>,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCNRQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,HAXX::quaternion<double>,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCTRQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,HAXX::quaternion<double>,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCCRQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,double,HAXX::quaternion<double>,'C','C'>();
}


BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNNCQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,HAXX::quaternion<double>,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNTCQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,HAXX::quaternion<double>,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNCCQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,HAXX::quaternion<double>,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTNCQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,HAXX::quaternion<double>,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTTCQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,HAXX::quaternion<double>,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTCCQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,HAXX::quaternion<double>,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCNCQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,HAXX::quaternion<double>,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCTCQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,HAXX::quaternion<double>,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCCCQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,std::complex<double>,HAXX::quaternion<double>,'C','C'>();
}




BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNNQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNTQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHNCQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTNQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTTQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHTCQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCNQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCTQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HHCCQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,HAXX::quaternion<double>,'C','C'>();
}





// Real-Quaternion matrix-matrix multiplication





BOOST_AUTO_TEST_CASE(hblas3_gemm_square_RHNNRR_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,double,double,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_RHNTRR_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,double,double,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_RHNCRR_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,double,double,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_RHTNRR_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,double,double,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_RHTTRR_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,double,double,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_RHTCRR_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,double,double,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_RHCNRR_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,double,double,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_RHCTRR_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,double,double,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_RHCCRR_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,double,double,'C','C'>();
}

BOOST_AUTO_TEST_SUITE_END()
