# Determine SIMD instruction set if any


message( STATUS "Attempting to determine SIMD instructions")

# LINUX checks
if( CMAKE_SYSTEM_NAME MATCHES "Linux" )

  # Get proc info
  EXEC_PROGRAM( cat ARGS "/proc/cpuinfo" OUTPUT_VARIABLE CPUINFO )


  # Check AVX and AVX2
  STRING( REGEX REPLACE "^.*(avx).*$" "\\1"  AVX_THERE  ${CPUINFO} )
  STRING( REGEX REPLACE "^.*(avx2).*$" "\\1" AVX2_THERE ${CPUINFO} )

  STRING( COMPARE EQUAL "avx"  "${AVX_THERE}"  AVX_TRUE  )
  STRING( COMPARE EQUAL "avx2" "${AVX2_THERE}" AVX2_TRUE )


# NON-LINUX defaults to generic code
else()

 set( AVX_TRUE    false)
 set( AVX2_TRUE   false)
 set( AVX512_TRUE false)

endif()


set( AVX_FOUND    ${AVX_TRUE}    CACHE BOOL "AVX instructions available"    )
set( AVX2_FOUND   ${AVX2_TRUE}   CACHE BOOL "AVX2 instructions available"   )
set( AVX512_FOUND ${AVX512_TRUE} CACHE BOOL "AVX512 instructions available" )


if( AVX2_FOUND )
  message( STATUS "-- AVX2 is largest available SIMD instruction set" )
elseif( AVX_FOUND )
  message( STATUS "-- AVX is largest available SIMD instruction set" )
endif()

