
add_subdirectory(fortran)

if( NOT ENABLE_GENERIC_FORTRAN )
  add_subdirectory(cxx)
endif()


configure_file( hblas1.cxx.in ${PROJECT_BINARY_DIR}/src/hblas/hblas1.cxx )
configure_file( hblas2.cxx.in ${PROJECT_BINARY_DIR}/src/hblas/hblas2.cxx )
configure_file( hblas3.cxx.in ${PROJECT_BINARY_DIR}/src/hblas/hblas3.cxx )

set( HBLAS_IMPL_SRC 
  ${PROJECT_BINARY_DIR}/src/hblas/hblas1.cxx
  ${PROJECT_BINARY_DIR}/src/hblas/hblas2.cxx
  ${PROJECT_BINARY_DIR}/src/hblas/hblas3.cxx
  util.cxx
)

add_library( hblas STATIC ${HBLAS_IMPL_SRC} )

target_link_libraries(hblas PUBLIC hblas_fortran)

if( NOT ENABLE_GENERIC_FORTRAN )
  target_link_libraries(hblas PUBLIC hblas_cxx)
endif()

install(TARGETS hblas DESTINATION lib)
