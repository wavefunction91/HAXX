# GEMM functions
add_library( hblas_gemmdd OBJECT gemm.cxx )
add_library( hblas_gemmdz OBJECT gemm.cxx )
add_library( hblas_gemmzd OBJECT gemm.cxx )
add_library( hblas_gemmzz OBJECT gemm.cxx )

target_compile_definitions( hblas_gemmdd PRIVATE "ALPHAF=DOUBLE"   "BETAF=DOUBLE"   )
target_compile_definitions( hblas_gemmdz PRIVATE "ALPHAF=DOUBLE"   "BETAF=DCOMPLEX" )
target_compile_definitions( hblas_gemmzd PRIVATE "ALPHAF=DCOMPLEX" "BETAF=DOUBLE"   )
target_compile_definitions( hblas_gemmzz PRIVATE "ALPHAF=DCOMPLEX" "BETAF=DCOMPLEX" )

set(GEMM_OBJ
    $<TARGET_OBJECTS:hblas_gemmdd>
    $<TARGET_OBJECTS:hblas_gemmzd>
    $<TARGET_OBJECTS:hblas_gemmdz>
    $<TARGET_OBJECTS:hblas_gemmzz>
  CACHE LIST "Object Files for GEMM functions"
)

