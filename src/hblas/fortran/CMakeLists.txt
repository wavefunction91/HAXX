
if( ENABLE_GENERIC_FORTRAN )
  set(HBLAS1_SRC 
    hblas1/hdotu.f hblas1/hdotc.f 
    hblas1/haxpydh.f hblas1/haxpyhh.f hblas1/haxpych.f 
    hblas1/hscalh.f hblas1/hscalc.f hblas1/hscald.f)
endif()

set(HBLAS2_SRC 
  hblas2/hgemvdd.f hblas2/hgemvzd.f hblas2/hgemvhd.f
  hblas2/hgemvdz.f hblas2/hgemvzz.f hblas2/hgemvhz.f
  hblas2/hgemvdh.f hblas2/hgemvzh.f hblas2/hgemvhh.f
  hblas2/hgerud.f hblas2/hgeruh.f hblas2/hgeruz.f
  hblas2/hgercd.f hblas2/hgerch.f hblas2/hgercz.f)
set(HBLAS3_SRC 
  hblas3/hgemmdd.f hblas3/hgemmzd.f hblas3/hgemmhd.f
  hblas3/hgemmdz.f hblas3/hgemmzz.f hblas3/hgemmhz.f
  hblas3/hgemmdh.f hblas3/hgemmzh.f hblas3/hgemmhh.f)

set(UTIL_SRC util/lsame.f util/hzexp.f util/hdexp.f util/hzcon.f)

add_library(hblas_fortran STATIC ${HBLAS1_SRC} ${HBLAS2_SRC} 
  ${HBLAS3_SRC} ${UTIL_SRC})


# Install FORTRAN library
install(TARGETS hblas_fortran DESTINATION lib)
