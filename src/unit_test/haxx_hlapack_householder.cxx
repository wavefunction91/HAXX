#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HLAPACK_HOUSEHOLDER

#include "haxx_ut.hpp"
#include "hblas/hblas1_impl.hpp"
#include "hblas/hblas2_impl.hpp"
#include "hblas/hblas3_impl.hpp"
#include "hlapack/householder_impl.hpp"


BOOST_AUTO_TEST_CASE(hlapack_larfg)
{
  std::vector<HAXX::quaternion<double>> 
    X(HBLAS1_VECLEN), XC(HBLAS1_VECLEN), SCR(HBLAS2_MATLEN), 
    SCR2(HBLAS2_MATLEN);

  
  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));


  std::copy(X.begin(),X.end(),XC.begin());

  double BETA = HLAPACK_LARFG(HBLAS1_VECLEN,&X[0],1);
  
  // Construct householder matrix
  std::fill(SCR.begin(),SCR.end(),0.);
  for(auto k = 0; k < HBLAS1_VECLEN; ++k) 
    SCR[RANK2_INDX(k,k,HBLAS1_VECLEN)] = 1.;

  HBLAS_GERC(HBLAS1_VECLEN,HBLAS1_VECLEN,-BETA,&X[0],1,&X[0],1,
    &SCR[0],HBLAS1_VECLEN);

  // Check Hermetian
  for(auto j = 0; j < HBLAS1_VECLEN; ++j)
  for(auto i = j; i < HBLAS1_VECLEN; ++i)
    BOOST_CHECK(
      CMP_Q(SCR[RANK2_INDX(i,j,HBLAS1_VECLEN)],
            HAXX::conj(SCR[RANK2_INDX(j,i,HBLAS1_VECLEN)]))
    );

  // Check Idempotent
  HBLAS_GEMM('N','N',HBLAS1_VECLEN,HBLAS1_VECLEN,HBLAS1_VECLEN,1.,
    &SCR[0],HBLAS1_VECLEN,&SCR[0],HBLAS1_VECLEN,0.,&SCR2[0],
    HBLAS1_VECLEN);

  for(auto j = 0; j < HBLAS1_VECLEN; ++j) {
    BOOST_CHECK(
      CMP_Q(SCR2[RANK2_INDX(j,j,HBLAS1_VECLEN)],HAXX::quaternion<double>(1.))
    );

    for(auto i = j+1; i < HBLAS1_VECLEN; ++i) {
      BOOST_CHECK_SMALL(SCR2[RANK2_INDX(i,j,HBLAS1_VECLEN)].real(),COMPARE_TOL);
      BOOST_CHECK_SMALL(SCR2[RANK2_INDX(i,j,HBLAS1_VECLEN)].imag_i(),
        COMPARE_TOL);
      BOOST_CHECK_SMALL(SCR2[RANK2_INDX(i,j,HBLAS1_VECLEN)].imag_j(),
        COMPARE_TOL);
      BOOST_CHECK_SMALL(SCR2[RANK2_INDX(i,j,HBLAS1_VECLEN)].imag_k(),
        COMPARE_TOL);

      BOOST_CHECK_SMALL(SCR2[RANK2_INDX(j,i,HBLAS1_VECLEN)].real(),COMPARE_TOL);
      BOOST_CHECK_SMALL(SCR2[RANK2_INDX(j,i,HBLAS1_VECLEN)].imag_i(),
        COMPARE_TOL);
      BOOST_CHECK_SMALL(SCR2[RANK2_INDX(j,i,HBLAS1_VECLEN)].imag_j(),
        COMPARE_TOL);
      BOOST_CHECK_SMALL(SCR2[RANK2_INDX(j,i,HBLAS1_VECLEN)].imag_k(),
        COMPARE_TOL);
    }
  }


  // Check for proper action on X
  auto Inner = HBLAS_DOTC(HBLAS1_VECLEN,&X[0],1,&XC[0],1);
  double x2nrm = HBLAS_NRM2(HBLAS1_VECLEN,&XC[0],1);
  auto   Uy    = XC[0] / HAXX::norm(XC[0]);

  HBLAS_SCAL('R',HBLAS1_VECLEN,Inner,&X[0],1);
  HBLAS_AXPY('L',HBLAS1_VECLEN,-BETA,&X[0],1,&XC[0],1);


  BOOST_CHECK(CMP_Q(-x2nrm*Uy,XC[0]));
  for(auto j = 1; j < HBLAS1_VECLEN; ++j){
    BOOST_CHECK_SMALL(XC[j].real()  ,COMPARE_TOL);
    BOOST_CHECK_SMALL(XC[j].imag_i(),COMPARE_TOL);
    BOOST_CHECK_SMALL(XC[j].imag_j(),COMPARE_TOL);
    BOOST_CHECK_SMALL(XC[j].imag_k(),COMPARE_TOL);
  }

  

}

BOOST_AUTO_TEST_CASE(hlapack_larf)
{
  std::vector<HAXX::quaternion<double>> 
    X(HBLAS1_VECLEN), XC(HBLAS1_VECLEN), SCR(HBLAS1_VECLEN);

  
  for(auto &x : X) 
    x = HAXX::quaternion<double>(dis(gen),dis(gen),dis(gen),dis(gen));

  std::copy(X.begin(),X.end(),XC.begin());

  HAXX::quaternion<double> BETA = HLAPACK_LARFG(HBLAS1_VECLEN,&X[0],1);
  auto Inner = HBLAS_DOTC(HBLAS1_VECLEN,&X[0],1,&XC[0],1);
  double x2nrm = HBLAS_NRM2(HBLAS1_VECLEN,&XC[0],1);
  auto   Uy    = XC[0] / HAXX::norm(XC[0]);

  HLAPACK_LARF('L',HBLAS1_VECLEN,1,&X[0],1,BETA,&XC[0],HBLAS1_VECLEN,
    &SCR[0]);

  BOOST_CHECK(CMP_Q(-x2nrm*Uy,XC[0]));
  for(auto j = 1; j < HBLAS1_VECLEN; ++j){
    BOOST_CHECK_SMALL(XC[j].real()  ,COMPARE_TOL);
    BOOST_CHECK_SMALL(XC[j].imag_i(),COMPARE_TOL);
    BOOST_CHECK_SMALL(XC[j].imag_j(),COMPARE_TOL);
    BOOST_CHECK_SMALL(XC[j].imag_k(),COMPARE_TOL);
  }
}
