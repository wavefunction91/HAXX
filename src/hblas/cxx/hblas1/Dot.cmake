# DOT functions

add_library( hblas_dotuv OBJECT dot.cxx )
add_library( hblas_dotcv OBJECT dot.cxx )

target_compile_definitions( hblas_dotuv PRIVATE "FNAME=HBLAS_DOTUV"         )
target_compile_definitions( hblas_dotcv PRIVATE "FNAME=HBLAS_DOTCV" "_CONJ" )

set( DOT_OBJ 
    $<TARGET_OBJECTS:hblas_dotuv> 
    $<TARGET_OBJECTS:hblas_dotcv> 
  CACHE LIST "Object Files for DOTV functions"
)

