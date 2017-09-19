

##### HBLAS1 COMPILATION #####

# DOT functions
configure_file(hblas1/dot.cxx.in 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotu.cxx)
configure_file(hblas1/dot.cxx.in 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotc.cxx)

set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotu.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_DOTU")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotc.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_DOTC" "_CONJ")


set( DOT_SRC 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotu.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotc.cxx 
)


# AXPY/SCAL functions
configure_file(hblas1/axpy_scal.cxx.in 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscald.cxx)
configure_file(hblas1/axpy_scal.cxx.in 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalz.cxx)
configure_file(hblas1/axpy_scal.cxx.in 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalh.cxx)

configure_file(hblas1/axpy_scal.cxx.in 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpydh.cxx)
configure_file(hblas1/axpy_scal.cxx.in 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyzh.cxx)
configure_file(hblas1/axpy_scal.cxx.in 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyhh.cxx)

set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscald.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_SCAL" "_SCAL" "ALPHAF=DOUBLE")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalz.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_SCAL" "_SCAL" "ALPHAF=DCOMPLEX")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalh.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_SCAL" "_SCAL" "ALPHAF=DQUATERNION")

set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpydh.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_AXPY" "_AXPY" "ALPHAF=DOUBLE")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyzh.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_AXPY" "_AXPY" "ALPHAF=DCOMPLEX")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyhh.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_AXPY" "_AXPY" "ALPHAF=DQUATERNION")


set(AXPY_SCAL_SRC
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscald.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalz.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalh.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpydh.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyzh.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyhh.cxx 
)

set(HBLAS1_SRC ${DOT_SRC} ${AXPY_SCAL_SRC})


##### HBLAS3 COMPILATION #####

# GEMM functions
configure_file(hblas3/gemm.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmdd.cxx)


set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmdd.cxx
  PROPERTY COMPILE_DEFINITIONS "ALPHAF=DOUBLE" "BETAF=DOUBLE")


set(GEMM_SRC
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmdd.cxx
)


set(HBLAS3_SRC ${GEMM_SRC})


##### HBLAS LIB TARGET #####

add_library(hblas_cxx STATIC ${HBLAS1_SRC} ${HBLAS3_SRC})




