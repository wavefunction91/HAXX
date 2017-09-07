/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HAXX_BLAS2

#include "haxx_ut.hpp"
#include "hblas/hblas1_impl.hpp"
#include "hblas/hblas2_def.hpp"


template <typename _AlphaF, typename _BetaF, char _TA>
void hblas2_gemv_square_LDSame() {

  // Random Quaternion vectors and matricies
  std::vector<HAXX::quaternion<double>> 
    A(HBLAS2_MATLEN), X(HBLAS1_VECLEN), Y(HBLAS1_VECLEN);

  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : Y) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));
  for(auto &x : A) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::vector<HAXX::quaternion<double>> YC(Y);

  _AlphaF ALPHA = genRandom<_AlphaF>();
  _BetaF  BETA = genRandom<_BetaF>();

  std::stringstream ss;

  ss << "hblas2_gemv_square_" << _TA;
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

  HBLAS_GEMV(_TA,HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&A[0],HBLAS1_VECLEN,&X[0],1,
    BETA,&Y[0],1);

  for(int i = 0; i < HBLAS1_VECLEN; i++) {
    HAXX_INT AROW = 0;
    HAXX_INT ACOL = i;
    HAXX_INT AStride = 1;
    if( _TA == 'N' ) {
      AROW = i;
      ACOL = 0;
      AStride = HBLAS1_VECLEN;
    }

    HAXX::quaternion<double> tmp;
    if( _TA == 'C' )
      tmp = HBLAS_DOTC(HBLAS1_VECLEN,&A[RANK2_INDX(AROW,ACOL,HBLAS1_VECLEN)],
          AStride,&X[0],1);
    else
      tmp = HBLAS_DOTU(HBLAS1_VECLEN,&A[RANK2_INDX(AROW,ACOL,HBLAS1_VECLEN)],
          AStride,&X[0],1);

    BOOST_CHECK(CMP_Q(ALPHA*tmp + BETA*YC[i],Y[i]));
  }

}

BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NRR_LDSame) {
  hblas2_gemv_square_LDSame<double,double,'N'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TRR_LDSame) {
  hblas2_gemv_square_LDSame<double,double,'T'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CRR_LDSame) {
  hblas2_gemv_square_LDSame<double,double,'C'>(); 
}


BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NCR_LDSame) {
  hblas2_gemv_square_LDSame<std::complex<double>,double,'N'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TCR_LDSame) {
  hblas2_gemv_square_LDSame<std::complex<double>,double,'T'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CCR_LDSame) {
  hblas2_gemv_square_LDSame<std::complex<double>,double,'C'>(); 
}


BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NQR_LDSame) {
  hblas2_gemv_square_LDSame<HAXX::quaternion<double>,double,'N'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TQR_LDSame) {
  hblas2_gemv_square_LDSame<HAXX::quaternion<double>,double,'T'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CQR_LDSame) {
  hblas2_gemv_square_LDSame<HAXX::quaternion<double>,double,'C'>(); 
}


BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NRC_LDSame) {
  hblas2_gemv_square_LDSame<double,std::complex<double>,'N'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TRC_LDSame) {
  hblas2_gemv_square_LDSame<double,std::complex<double>,'T'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CRC_LDSame) {
  hblas2_gemv_square_LDSame<double,std::complex<double>,'C'>(); 
}


BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NCC_LDSame) {
  hblas2_gemv_square_LDSame<std::complex<double>,std::complex<double>,'N'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TCC_LDSame) {
  hblas2_gemv_square_LDSame<std::complex<double>,std::complex<double>,'T'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CCC_LDSame) {
  hblas2_gemv_square_LDSame<std::complex<double>,std::complex<double>,'C'>(); 
}


BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NQC_LDSame) {
  hblas2_gemv_square_LDSame<HAXX::quaternion<double>,std::complex<double>,'N'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TQC_LDSame) {
  hblas2_gemv_square_LDSame<HAXX::quaternion<double>,std::complex<double>,'T'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CQC_LDSame) {
  hblas2_gemv_square_LDSame<HAXX::quaternion<double>,std::complex<double>,'C'>(); 
}



BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NRQ_LDSame) {
  hblas2_gemv_square_LDSame<double,HAXX::quaternion<double>,'N'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TRQ_LDSame) {
  hblas2_gemv_square_LDSame<double,HAXX::quaternion<double>,'T'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CRQ_LDSame) {
  hblas2_gemv_square_LDSame<double,HAXX::quaternion<double>,'C'>(); 
}


BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NCQ_LDSame) {
  hblas2_gemv_square_LDSame<std::complex<double>,HAXX::quaternion<double>,'N'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TCQ_LDSame) {
  hblas2_gemv_square_LDSame<std::complex<double>,HAXX::quaternion<double>,'T'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CCQ_LDSame) {
  hblas2_gemv_square_LDSame<std::complex<double>,HAXX::quaternion<double>,'C'>(); 
}


BOOST_AUTO_TEST_CASE(hblas2_gemv_square_NQQ_LDSame) {
  hblas2_gemv_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,'N'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_TQQ_LDSame) {
  hblas2_gemv_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,'T'>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gemv_square_CQQ_LDSame) {
  hblas2_gemv_square_LDSame<HAXX::quaternion<double>,HAXX::quaternion<double>,'C'>(); 
}





























template <typename _AlphaF, bool _CONJ>
void hblas2_ger_square() {

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

  _AlphaF ALPHA = genRandom<_AlphaF>();

  std::stringstream ss;

  ss << "hblas2_ger ";
  if( _CONJ ) ss << "c";
  else        ss << "u";
  ss << "_square_";
  if(std::is_same<double,_AlphaF>::value) ss << "R";
  else if(std::is_same<std::complex<double>,_AlphaF>::value) ss << "C";
  else ss << "Q";

  if( _CONJ )
    HBLAS_GERC(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);
  else
    HBLAS_GERU(HBLAS1_VECLEN,HBLAS1_VECLEN,ALPHA,&X[0],1,&Y[0],1,&A[0],HBLAS1_VECLEN);

  for(int j = 0; j < HBLAS1_VECLEN; j++) 
  for(int i = 0; i < HBLAS1_VECLEN; i++) {

    if( _CONJ )
      BOOST_CHECK(
        CMP_Q(
          ALPHA*X[i]*HAXX::conj(Y[j]) + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          A[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );
    else
      BOOST_CHECK(
        CMP_Q(
          ALPHA*X[i]*Y[j] + AC[RANK2_INDX(i,j,HBLAS1_VECLEN)],
          A[RANK2_INDX(i,j,HBLAS1_VECLEN)]
        )
      );

  }
}




BOOST_AUTO_TEST_CASE(hblas2_geru_square_R) {
  hblas2_ger_square<double,false>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_geru_square_C) {
  hblas2_ger_square<std::complex<double>,false>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_geru_square_Q) {
  hblas2_ger_square<HAXX::quaternion<double>,false>(); 
}

BOOST_AUTO_TEST_CASE(hblas2_gerc_square_R) {
  hblas2_ger_square<double,true>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gerc_square_C) {
  hblas2_ger_square<std::complex<double>,true>(); 
}
BOOST_AUTO_TEST_CASE(hblas2_gerc_square_Q) {
  hblas2_ger_square<HAXX::quaternion<double>,true>(); 
}
