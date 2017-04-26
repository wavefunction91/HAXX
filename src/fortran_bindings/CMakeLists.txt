set(HBLAS1_SRC hblas/hblas1/hdotu.f hblas/hblas1/hdotc.f hblas/hblas1/haxpydh.f
  hblas/hblas1/haxpyhh.f hblas/hblas1/haxpych.f hblas/hblas1/hscalh.f
  hblas/hblas1/hscalc.f hblas/hblas1/hscald.f)
set(HBLAS2_SRC hblas/hblas2/hgemvdd.f)
set(HBLAS3_SRC hblas/hblas3/hgemmdd.f)

add_library(haxx_fortran STATIC ${HBLAS1_SRC} ${HBLAS2_SRC} 
  ${HBLAS3_SRC} util/lsame.f)
