##### HBLAS3 COMPILATION #####

include(Gemm.cmake)

set(HBLAS3_OBJ ${GEMM_OBJ} CACHE LIST 
  "Object Files for HBLAS3 functions" )
