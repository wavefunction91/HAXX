
# COPY / SWAP functions

configure_file(hblas1/copy_swap.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hcopy.cxx)
configure_file(hblas1/copy_swap.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hswap.cxx)

set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hcopy.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_COPY" "_COPY")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hswap.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_SWAP" "_SWAP")


set( COPY_SWAP_SRC 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hcopy.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hswap.cxx 
)

# DOT functions
configure_file(hblas1/dot.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotu.cxx)
configure_file(hblas1/dot.cxx 
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
configure_file(hblas1/axpy_scal.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscald.cxx)
configure_file(hblas1/axpy_scal.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalz.cxx)
configure_file(hblas1/axpy_scal.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalh.cxx)

configure_file(hblas1/axpy_scal.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpydh.cxx)
configure_file(hblas1/axpy_scal.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyzh.cxx)
configure_file(hblas1/axpy_scal.cxx 
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

set(HBLAS1_SRC ${DOT_SRC} ${AXPY_SCAL_SRC} ${COPY_SWAP_SRC})

add_library(hblas_cxx STATIC ${HBLAS1_SRC})
