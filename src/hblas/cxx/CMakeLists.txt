add_subdirectory(hblas1)       
add_subdirectory(hblas3)       

##### HBLAS LIB TARGET #####

add_library(hblas_cxx STATIC ${HBLAS1_OBJ} ${HBLAS3_OBJ})
target_include_directories(hblas_cxx PUBLIC ${PROJECT_SOURCE_DIR}/src/hblas/cxx/hblas3)




