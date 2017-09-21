

##### HBLAS1 COMPILATION #####
       
# COPY / SWAP functions

configure_file(hblas1/copy_swap.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hcopyv.cxx)
configure_file(hblas1/copy_swap.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hswapv.cxx)

set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hcopyv.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_COPYV" "_COPY")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hswapv.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_SWAPV" "_SWAP")


set( COPY_SWAP_SRC 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hcopyv.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hswapv.cxx 
)
              

# DOT functions
configure_file(hblas1/dot.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotuv.cxx)
configure_file(hblas1/dot.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotcv.cxx)

set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotuv.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_DOTUV")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotcv.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_DOTCV" "_CONJ")


set( DOT_SRC 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotuv.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hdotcv.cxx 
)


# AXPY/SCAL functions
configure_file(hblas1/axpy_scal.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalvd.cxx)
configure_file(hblas1/axpy_scal.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalvz.cxx)
configure_file(hblas1/axpy_scal.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalvh.cxx)

configure_file(hblas1/axpy_scal.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyvdh.cxx)
configure_file(hblas1/axpy_scal.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyvzh.cxx)
configure_file(hblas1/axpy_scal.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyvhh.cxx)

set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalvd.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_SCALV" "_SCAL" "ALPHAF=DOUBLE")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalvz.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_SCALV" "_SCAL" "ALPHAF=DCOMPLEX")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalvh.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_SCALV" "_SCAL" "ALPHAF=DQUATERNION")

set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyvdh.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_AXPYV" "_AXPY" "ALPHAF=DOUBLE")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyvzh.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_AXPYV" "_AXPY" "ALPHAF=DCOMPLEX")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyvhh.cxx 
  PROPERTY COMPILE_DEFINITIONS "FNAME=HBLAS_AXPYV" "_AXPY" "ALPHAF=DQUATERNION")


set(AXPY_SCAL_SRC
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalvd.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalvz.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/hscalvh.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyvdh.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyvzh.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas1/haxpyvhh.cxx 
)

set(HBLAS1_SRC ${DOT_SRC} ${AXPY_SCAL_SRC} ${COPY_SWAP_SRC})


##### HBLAS3 COMPILATION #####

# GEMM functions
configure_file(hblas3/gemm.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmdd.cxx)
configure_file(hblas3/gemm.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmdz.cxx)
configure_file(hblas3/gemm.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmzd.cxx)
configure_file(hblas3/gemm.cxx 
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmzz.cxx)


set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmdd.cxx
  PROPERTY COMPILE_DEFINITIONS "ALPHAF=DOUBLE" "BETAF=DOUBLE")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmdz.cxx
  PROPERTY COMPILE_DEFINITIONS "ALPHAF=DOUBLE" "BETAF=DCOMPLEX")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmzd.cxx
  PROPERTY COMPILE_DEFINITIONS "ALPHAF=DCOMPLEX" "BETAF=DOUBLE")
set_property(SOURCE ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmzz.cxx
  PROPERTY COMPILE_DEFINITIONS "ALPHAF=DCOMPLEX" "BETAF=DCOMPLEX")


set(GEMM_SRC
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmdd.cxx
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmzd.cxx
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmdz.cxx
  ${PROJECT_BINARY_DIR}/src/hblas/cxx/hblas3/hgemmzz.cxx
)

set(HBLAS3_SRC ${GEMM_SRC})


##### HBLAS LIB TARGET #####

add_library(hblas_cxx STATIC ${HBLAS1_SRC} ${HBLAS3_SRC})
target_include_directories(hblas_cxx PUBLIC ${PROJECT_SOURCE_DIR}/src/hblas/cxx/hblas3)




