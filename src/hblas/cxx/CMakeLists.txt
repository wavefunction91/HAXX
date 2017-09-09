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

set(HBLAS1_SRC ${DOT_SRC})

add_library(hblas_cxx STATIC ${HBLAS1_SRC})
