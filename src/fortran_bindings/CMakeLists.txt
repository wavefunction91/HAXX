set(HBLAS1_SRC hblas/hblas1/hdotu.f hblas/hblas1/hdotc.f hblas/hblas1/haxpydh.f
  hblas/hblas1/haxpyhh.f hblas/hblas1/haxpych.f hblas/hblas1/hscalh.f
  hblas/hblas1/hscalc.f hblas/hblas1/hscald.f)
set(HBLAS2_SRC 
  hblas/hblas2/hgemvdd.f hblas/hblas2/hgemvzd.f hblas/hblas2/hgemvhd.f
  hblas/hblas2/hgemvdz.f hblas/hblas2/hgemvzz.f hblas/hblas2/hgemvhz.f
  hblas/hblas2/hgemvdh.f hblas/hblas2/hgemvzh.f hblas/hblas2/hgemvhh.f
  hblas/hblas2/hgerud.f hblas/hblas2/hgeruh.f hblas/hblas2/hgeruz.f
  hblas/hblas2/hgercd.f hblas/hblas2/hgerch.f hblas/hblas2/hgercz.f)
set(HBLAS3_SRC 
  hblas/hblas3/hgemmdd.f hblas/hblas3/hgemmzd.f hblas/hblas3/hgemmhd.f
  hblas/hblas3/hgemmdz.f hblas/hblas3/hgemmzz.f hblas/hblas3/hgemmhz.f
  hblas/hblas3/hgemmdh.f hblas/hblas3/hgemmzh.f hblas/hblas3/hgemmhh.f)

set(UTIL_SRC util/lsame.f util/hzexp.f util/hdexp.f util/hzcon.f)

add_library(haxx_fortran STATIC ${HBLAS1_SRC} ${HBLAS2_SRC} 
  ${HBLAS3_SRC} ${UTIL_SRC})


# Install FORTRAN library
install(TARGETS haxx_fortran DESTINATION lib)
