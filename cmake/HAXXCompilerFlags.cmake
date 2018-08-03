include(CheckCXXCompilerFlag)
include(CheckFortranCompilerFlag)

# Handle C++14 Flags
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
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


# Check IPO
check_cxx_compiler_flag("-ipo" CXX_USES_IPO)

if( CXX_USES_IPO AND NOT DISABLE_IPO )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ipo" )
endif()


# Host SIMD flags
if( HAXX_USE_HOST_SIMD )

  # Determine CXX opt flags
  check_cxx_compiler_flag("-march=native" CXX_USES_MARCH_NATIVE)
  check_cxx_compiler_flag("-xHost"        CXX_USES_XHOST       )

  check_fortran_compiler_flag("-march=native" FC_USES_MARCH_NATIVE)
  check_fortran_compiler_flag("-xHost"        FC_USES_XHOST       )

  # Add Host flags
  if( CXX_USES_XHOST )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost" )
  elseif( CXX_USES_MARCH_NATIVE )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3" )
  else()
    message( WARNING "Unable to determine proper HOST flags for CXX compiler" )
  endif()

  if( FC_USES_XHOST )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost" )
  elseif( FC_USES_MARCH_NATIVE )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native" )
  else()
    message( WARNING "Unable to determine proper HOST flags for FC" )
  endif()

endif()

# HAXX Types

# Index Integer Type
if( NOT HAXX_INT )
  set( HAXX_INT int32_t )
endif()
add_definitions("-DHAXX_INT=${HAXX_INT}")


