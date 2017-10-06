# AXPY/SCAL functions

add_library( hblas_scalvd  OBJECT axpy_scal.cxx )
add_library( hblas_scalvz  OBJECT axpy_scal.cxx )
add_library( hblas_scalvh  OBJECT axpy_scal.cxx )
add_library( hblas_axpyvdh OBJECT axpy_scal.cxx )
add_library( hblas_axpyvzh OBJECT axpy_scal.cxx )
add_library( hblas_axpyvhh OBJECT axpy_scal.cxx )

target_compile_definitions( hblas_scalvd  PRIVATE "FNAME=HBLAS_SCALV" "_SCAL" "ALPHAF=DOUBLE"      )
target_compile_definitions( hblas_scalvz  PRIVATE "FNAME=HBLAS_SCALV" "_SCAL" "ALPHAF=DCOMPLEX"    )
target_compile_definitions( hblas_scalvh  PRIVATE "FNAME=HBLAS_SCALV" "_SCAL" "ALPHAF=DQUATERNION" )
target_compile_definitions( hblas_axpyvdh PRIVATE "FNAME=HBLAS_AXPYV" "_AXPY" "ALPHAF=DOUBLE"      )
target_compile_definitions( hblas_axpyvzh PRIVATE "FNAME=HBLAS_AXPYV" "_AXPY" "ALPHAF=DCOMPLEX"    )
target_compile_definitions( hblas_axpyvhh PRIVATE "FNAME=HBLAS_AXPYV" "_AXPY" "ALPHAF=DQUATERNION" )

set(AXPY_SCAL_OBJ
    $<TARGET_OBJECTS:hblas_scalvd> 
    $<TARGET_OBJECTS:hblas_scalvz> 
    $<TARGET_OBJECTS:hblas_scalvh> 
    $<TARGET_OBJECTS:hblas_axpyvdh> 
    $<TARGET_OBJECTS:hblas_axpyvzh> 
    $<TARGET_OBJECTS:hblas_axpyvhh> 
  CACHE LIST "Object Files for SCAL and AXPY functions"
)

