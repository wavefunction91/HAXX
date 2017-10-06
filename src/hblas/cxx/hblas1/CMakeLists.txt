##### HBLAS1 COMPILATION #####

include(CopySwap.cmake)
include(Dot.cmake)
include(AXPYScal.cmake)


set( HBLAS1_OBJ ${DOT_OBJ} ${AXPY_SCAL_OBJ} ${COPY_SWAP_OBJ}
  CACHE LIST "Object Files for HBLAS1 functions" )
