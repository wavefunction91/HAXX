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

template <typename _AlphaF, typename _BetaF, char _TA, char _TB>
void hblas3_gemm_square_LDSame() {

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

  _AlphaF ALPHA = genRandom<_AlphaF>();
  _BetaF  BETA = genRandom<_BetaF>();

  std::stringstream ss;

  ss << "hblas3_gemm_square_" << _TA << _TB;
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
      BOOST_CHECK( 
        CMP_Q(
          SCR[i] + BETA*CC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          C[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    }
  }
};


BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNRR_LDSame) {
  hblas3_gemm_square_LDSame<double,double,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTRR_LDSame) {
  hblas3_gemm_square_LDSame<double,double,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCRR_LDSame) {
  hblas3_gemm_square_LDSame<double,double,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNRR_LDSame) {
  hblas3_gemm_square_LDSame<double,double,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTRR_LDSame) {
  hblas3_gemm_square_LDSame<double,double,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCRR_LDSame) {
  hblas3_gemm_square_LDSame<double,double,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNRR_LDSame) {
  hblas3_gemm_square_LDSame<double,double,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTRR_LDSame) {
  hblas3_gemm_square_LDSame<double,double,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCRR_LDSame) {
  hblas3_gemm_square_LDSame<double,double,'C','C'>();
}


BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNCR_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,double,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTCR_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,double,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCCR_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,double,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNCR_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,double,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTCR_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,double,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCCR_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,double,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNCR_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,double,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTCR_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,double,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCCR_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,double,'C','C'>();
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,double,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,double,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,double,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,double,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,double,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,double,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,double,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,double,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCQR_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,double,'C','C'>();
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNRC_LDSame) {
  hblas3_gemm_square_LDSame<double,std::complex<double>,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTRC_LDSame) {
  hblas3_gemm_square_LDSame<double,std::complex<double>,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCRC_LDSame) {
  hblas3_gemm_square_LDSame<double,std::complex<double>,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNRC_LDSame) {
  hblas3_gemm_square_LDSame<double,std::complex<double>,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTRC_LDSame) {
  hblas3_gemm_square_LDSame<double,std::complex<double>,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCRC_LDSame) {
  hblas3_gemm_square_LDSame<double,std::complex<double>,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNRC_LDSame) {
  hblas3_gemm_square_LDSame<double,std::complex<double>,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTRC_LDSame) {
  hblas3_gemm_square_LDSame<double,std::complex<double>,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCRC_LDSame) {
  hblas3_gemm_square_LDSame<double,std::complex<double>,'C','C'>();
}


BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNCC_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,std::complex<double>,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTCC_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,std::complex<double>,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCCC_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,std::complex<double>,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNCC_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,std::complex<double>,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTCC_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,std::complex<double>,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCCC_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,std::complex<double>,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNCC_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,std::complex<double>,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTCC_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,std::complex<double>,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCCC_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,std::complex<double>,'C','C'>();
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,std::complex<double>,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,std::complex<double>,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,std::complex<double>,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,std::complex<double>,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,std::complex<double>,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,std::complex<double>,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,std::complex<double>,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,std::complex<double>,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCQC_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,std::complex<double>,'C','C'>();
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNRQ_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTRQ_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCRQ_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNRQ_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTRQ_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCRQ_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNRQ_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTRQ_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCRQ_LDSame) {
  hblas3_gemm_square_LDSame<double,HAXX::quaternion<double>,'C','C'>();
}


BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNCQ_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,HAXX::quaternion<double>,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTCQ_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,HAXX::quaternion<double>,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCCQ_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,HAXX::quaternion<double>,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNCQ_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,HAXX::quaternion<double>,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTCQ_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,HAXX::quaternion<double>,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCCQ_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,HAXX::quaternion<double>,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNCQ_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,HAXX::quaternion<double>,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTCQ_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,HAXX::quaternion<double>,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCCQ_LDSame) {
  hblas3_gemm_square_LDSame<std::complex<double>,HAXX::quaternion<double>,'C','C'>();
}

BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NNQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NTQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_NCQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TNQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TTQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_TCQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CNQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CTQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_CCQQ_LDSame) {
  hblas3_gemm_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,'C','C'>();
}
