# COPY / SWAP functions

add_library( hblas_copyv OBJECT copy_swap.cxx )
add_library( hblas_swapv OBJECT copy_swap.cxx )
              
target_compile_definitions( hblas_copyv PRIVATE "FNAME=HBLAS_COPYV" "_COPY" )
target_compile_definitions( hblas_swapv PRIVATE "FNAME=HBLAS_SWAPV" "_SWAP" )

set( COPY_SWAP_OBJ 
    $<TARGET_OBJECTS:hblas_copyv> 
    $<TARGET_OBJECTS:hblas_swapv> 
  CACHE LIST "Object Files for COPYV and SWAPV function" )


