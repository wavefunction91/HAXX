include(CheckCXXCompilerFlag)
include(CheckFortranCompilerFlag)

# Handle C++11 Flags
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

# Check for FOTRAN preprocessor
check_fortran_compiler_flag("-fpp" FC_USES_FPP)
check_fortran_compiler_flag("-cpp" FC_USES_CPP)

if(FC_USES_FPP)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp")
elseif(FC_USES_CPP)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
else()
  message(FATAL "Unable to Determine a Suitable FORTRAN Preprocessor")
endif()
