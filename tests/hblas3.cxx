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


// Typedefs
using d  = double;
using cd = std::complex<d>;
using qd = HAXX::quaternion<d>;

template< class T, class U >
inline constexpr bool is_same_v = std::is_same<T, U>::value;

template< class T >
inline std::string t_str() {

  return
  (is_same_v<T,  d>) ? "DOUBLE" :
  (is_same_v<T, cd>) ? "DOUBLE COMPLEX" :
                       "DOUBLE QUATERNION";
}


BOOST_AUTO_TEST_SUITE(HBLAS3)




template <typename _AF, typename _BF, typename _AlphaF, typename _BetaF, 
  char _TA, char _TB, size_t N = HBLAS1_VECLEN, size_t M = HBLAS1_VECLEN, 
  size_t K = HBLAS1_VECLEN >
void hblas3_gemm() {

  size_t NK = N*K;
  size_t KM = K*M;
  size_t NM = N*M;

  size_t LDA = ( _TA == 'N' ) ? N : K;
  size_t LDB = ( _TB == 'N' ) ? K : M;
  size_t LDC = N;

  // Random Quaternion vectors and matricies
  std::vector<_AF> A(NK);
  std::vector<_BF> B(KM);
  std::vector<qd>  C(NM);

  for(auto &x : A) x = genRandom<_AF>();
  for(auto &x : B) x = genRandom<_BF>();
  for(auto &x : C) x = genRandom<qd> ();

  std::vector<qd> CC(C);
  std::vector<qd> SCR(N);

  _AlphaF ALPHA = genRandom<_AlphaF>();
  _BetaF  BETA  = genRandom<_BetaF>();


  std::string AFstr     = t_str<_AF>(); 
  std::string BFstr     = t_str<_BF>(); 
  std::string AlphaFstr = t_str<_AlphaF>(); 
  std::string BetaFstr  = t_str<_BetaF>(); 
                
                


  std::cout << "\n HBLAS GEMM TEST: \n";
  std::cout << "   - _F =      DOUBLE    \n";
  std::cout << "   - _AF =     " << AFstr << "\n";
  std::cout << "   - _BF =     " << BFstr << "\n";
  std::cout << "   - _AlphaF = " << AlphaFstr << "\n";
  std::cout << "   - _BetaF =  " << BetaFstr << "\n";
  std::cout << "   - _TA =     " << _TA << "\n";
  std::cout << "   - _TB =     " << _TB << "\n";

  std::cout << "\n   PROBLEM DIMENSIONS\n";
  std::cout << "     - N =     " << N << "\n";
  std::cout << "     - M =     " << M << "\n";
  std::cout << "     - K =     " << K << "\n";
  std::cout << "     - LDA =   " << LDA << "\n";
  std::cout << "     - LDB =   " << LDB << "\n";
  std::cout << "     - LDC =   " << LDC << "\n";



  std::cout << "\n";



  auto gemmStart = std::chrono::high_resolution_clock::now();

  HBLAS_GEMM(_TA, _TB, N, M, K, ALPHA, &A[0], LDA, &B[0], LDB, BETA, 
    &C[0], LDC);

  auto gemmEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> gemmDur = gemmEnd - gemmStart;
  std::cout << std::endl << " * GEMM took " << gemmDur.count() 
    << " s\n\n";


  HAXX_INT BStride = 1;
  if( _TB != 'N' ) BStride = LDB;
  if( _TB == 'C' ) 
    for(auto &x : B) x = HAXX::conj(x); 

  for(size_t j = 0; j < M; j++) {

    HAXX_INT BStart = j;
    if( _TB == 'N' ) BStart *= LDB;

    HBLAS_GEMV(_TA, N, K, ALPHA, &A[0], LDA, &B[BStart], BStride, 0., 
      &SCR[0], 1);

    for(size_t i = 0; i < N; i++) {

      size_t indx = RANK2_INDX(i,j,LDC);
      auto refVal = SCR[i] + BETA*CC[indx];
      auto tstVal = C[indx];

      BOOST_CHECK_MESSAGE( 
        CMP_Q( refVal, tstVal ),
        "\n i = " << i << 
          " j = " << j << 
          ", HBLAS_GEMM = " << tstVal << 
          ", REFERENCE = "  << refVal 
      );
    }
  }

};



// Exhaustive Tests...


// Quaternion-Quaternion matrix-matrix multiplication


// Double ALPHA / Double BETA
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NN_DD) {
  hblas3_gemm<qd,qd,d,d,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NT_DD) {
  hblas3_gemm<qd,qd,d,d,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NC_DD) {
  hblas3_gemm<qd,qd,d,d,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TN_DD) {
  hblas3_gemm<qd,qd,d,d,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TT_DD) {
  hblas3_gemm<qd,qd,d,d,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TC_DD) {
  hblas3_gemm<qd,qd,d,d,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CN_DD) {
  hblas3_gemm<qd,qd,d,d,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CT_DD) {
  hblas3_gemm<qd,qd,d,d,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CC_DD) {
  hblas3_gemm<qd,qd,d,d,'C','C'>();
}

// Double Complex ALPHA / Double BETA
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NN_ZD) {
  hblas3_gemm<qd,qd,cd,d,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NT_ZD) {
  hblas3_gemm<qd,qd,cd,d,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NC_ZD) {
  hblas3_gemm<qd,qd,cd,d,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TN_ZD) {
  hblas3_gemm<qd,qd,cd,d,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TT_ZD) {
  hblas3_gemm<qd,qd,cd,d,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TC_ZD) {
  hblas3_gemm<qd,qd,cd,d,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CN_ZD) {
  hblas3_gemm<qd,qd,cd,d,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CT_ZD) {
  hblas3_gemm<qd,qd,cd,d,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CC_ZD) {
  hblas3_gemm<qd,qd,cd,d,'C','C'>();
}



// Double Quaternion ALPHA / Double BETA
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NN_HD) {
  hblas3_gemm<qd,qd,qd,d,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NT_HD) {
  hblas3_gemm<qd,qd,qd,d,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NC_HD) {
  hblas3_gemm<qd,qd,qd,d,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TN_HD) {
  hblas3_gemm<qd,qd,qd,d,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TT_HD) {
  hblas3_gemm<qd,qd,qd,d,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TC_HD) {
  hblas3_gemm<qd,qd,qd,d,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CN_HD) {
  hblas3_gemm<qd,qd,qd,d,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CT_HD) {
  hblas3_gemm<qd,qd,qd,d,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CC_HD) {
  hblas3_gemm<qd,qd,qd,d,'C','C'>();
}

// Double ALPHA / Double Complex BETA
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NN_DZ) {
  hblas3_gemm<qd,qd,d,cd,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NT_DZ) {
  hblas3_gemm<qd,qd,d,cd,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NZ_DC) {
  hblas3_gemm<qd,qd,d,cd,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TN_DZ) {
  hblas3_gemm<qd,qd,d,cd,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TT_DZ) {
  hblas3_gemm<qd,qd,d,cd,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TZ_DC) {
  hblas3_gemm<qd,qd,d,cd,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CN_DZ) {
  hblas3_gemm<qd,qd,d,cd,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CT_DZ) {
  hblas3_gemm<qd,qd,d,cd,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CZ_DC) {
  hblas3_gemm<qd,qd,d,cd,'C','C'>();
}


// Double Complex ALPHA / Double Complex BETA
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NN_ZZ) {
  hblas3_gemm<qd,qd,cd,cd,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NT_ZZ) {
  hblas3_gemm<qd,qd,cd,cd,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NC_ZZ) {
  hblas3_gemm<qd,qd,cd,cd,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TN_ZZ) {
  hblas3_gemm<qd,qd,cd,cd,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TT_ZZ) {
  hblas3_gemm<qd,qd,cd,cd,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TC_ZZ) {
  hblas3_gemm<qd,qd,cd,cd,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CN_ZZ) {
  hblas3_gemm<qd,qd,cd,cd,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CT_ZZ) {
  hblas3_gemm<qd,qd,cd,cd,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CC_ZZ) {
  hblas3_gemm<qd,qd,cd,cd,'C','C'>();
}


// Double Quaternion ALPHA / Double Complex BETA
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NN_HZ) {
  hblas3_gemm<qd,qd,qd,cd,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NT_HZ) {
  hblas3_gemm<qd,qd,qd,cd,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NC_HZ) {
  hblas3_gemm<qd,qd,qd,cd,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TN_HZ) {
  hblas3_gemm<qd,qd,qd,cd,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TT_HZ) {
  hblas3_gemm<qd,qd,qd,cd,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TC_HZ) {
  hblas3_gemm<qd,qd,qd,cd,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CN_HZ) {
  hblas3_gemm<qd,qd,qd,cd,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CT_HZ) {
  hblas3_gemm<qd,qd,qd,cd,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CC_HZ) {
  hblas3_gemm<qd,qd,qd,cd,'C','C'>();
}

// Double ALPHA / Double Quaternion BETA
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NN_DH) {
  hblas3_gemm<qd,qd,d,qd,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NT_DH) {
  hblas3_gemm<qd,qd,d,qd,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NC_DH) {
  hblas3_gemm<qd,qd,d,qd,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TN_DH) {
  hblas3_gemm<qd,qd,d,qd,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TT_DH) {
  hblas3_gemm<qd,qd,d,qd,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TC_DH) {
  hblas3_gemm<qd,qd,d,qd,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CN_DH) {
  hblas3_gemm<qd,qd,d,qd,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CT_DH) {
  hblas3_gemm<qd,qd,d,qd,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CC_DH) {
  hblas3_gemm<qd,qd,d,qd,'C','C'>();
}


// Double Complex ALPHA / Double Quaternion BETA
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NN_ZH) {
  hblas3_gemm<qd,qd,cd,qd,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NT_ZH) {
  hblas3_gemm<qd,qd,cd,qd,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NC_ZH) {
  hblas3_gemm<qd,qd,cd,qd,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TN_ZH) {
  hblas3_gemm<qd,qd,cd,qd,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TT_ZH) {
  hblas3_gemm<qd,qd,cd,qd,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TC_ZH) {
  hblas3_gemm<qd,qd,cd,qd,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CN_ZH) {
  hblas3_gemm<qd,qd,cd,qd,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CT_ZH) {
  hblas3_gemm<qd,qd,cd,qd,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CC_ZH) {
  hblas3_gemm<qd,qd,cd,qd,'C','C'>();
}


// Double Quaternion ALPHA / Double Quaternion BETA
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NN_HH) {
  hblas3_gemm<qd,qd,qd,qd,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NT_HH) {
  hblas3_gemm<qd,qd,qd,qd,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_NC_HH) {
  hblas3_gemm<qd,qd,qd,qd,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TN_HH) {
  hblas3_gemm<qd,qd,qd,qd,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TT_HH) {
  hblas3_gemm<qd,qd,qd,qd,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_TC_HH) {
  hblas3_gemm<qd,qd,qd,qd,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CN_HH) {
  hblas3_gemm<qd,qd,qd,qd,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CT_HH) {
  hblas3_gemm<qd,qd,qd,qd,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_HH_CC_HH) {
  hblas3_gemm<qd,qd,qd,qd,'C','C'>();
}





// Real-Quaternion matrix-matrix multiplication

// Double ALPHA / Double BETA
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_DH_NN_DD) {
  hblas3_gemm<d,qd,d,d,'N','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_DH_NT_DD) {
  hblas3_gemm<d,qd,d,d,'N','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_DH_NC_DD) {
  hblas3_gemm<d,qd,d,d,'N','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_DH_TN_DD) {
  hblas3_gemm<d,qd,d,d,'T','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_DH_TT_DD) {
  hblas3_gemm<d,qd,d,d,'T','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_DH_TC_DD) {
  hblas3_gemm<d,qd,d,d,'T','C'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_DH_CN_DD) {
  hblas3_gemm<d,qd,d,d,'C','N'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_DH_CT_DD) {
  hblas3_gemm<d,qd,d,d,'C','T'>();
}
BOOST_AUTO_TEST_CASE(hblas3_gemm_square_DH_CC_DD) {
  hblas3_gemm<d,qd,d,d,'C','C'>();
}


// Misc fringe cases
//

BOOST_AUTO_TEST_CASE(hblas_gemm_N_NSQ) {
  hblas3_gemm<qd,qd,d,d,'N','N',66>();
};
BOOST_AUTO_TEST_CASE(hblas_gemm_M_NSQ) {
  hblas3_gemm<qd,qd,d,d,'N','N',HBLAS1_VECLEN,66>();
};
BOOST_AUTO_TEST_CASE(hblas_gemm_K_NSQ) {
  hblas3_gemm<qd,qd,d,d,'N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,66>();
};


BOOST_AUTO_TEST_CASE(hblas_gemm_N_ODD) {
  hblas3_gemm<qd,qd,d,d,'N','N',67>();
};
BOOST_AUTO_TEST_CASE(hblas_gemm_M_ODD) {
  hblas3_gemm<qd,qd,d,d,'N','N',HBLAS1_VECLEN,67>();
};
BOOST_AUTO_TEST_CASE(hblas_gemm_K_ODD) {
  hblas3_gemm<qd,qd,d,d,'N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,67>();
};

BOOST_AUTO_TEST_SUITE_END()
