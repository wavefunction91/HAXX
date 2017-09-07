
add_subdirectory(fortran_bindings)

add_library(hblas STATIC hblas3.cxx)
target_link_libraries(hblas PUBLIC hblas_fortran)

install(TARGETS hblas DESTINATION lib)
