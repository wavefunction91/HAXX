if( HAXX_USE_HOST_SIMD )

  # Try to determine SIMD 
  include( FindSIMD )
  
  if( NOT AVX_FOUND AND NOT AVX2_FOUND AND NOT AVX512_FOUND )
    message( WARNING "HAXX only provided optimal implementations for AVX, AVX2 and AVX-512 -- Defaulting to Generic FORTRAN build" )
  
    set( ENABLE_GENERIC_FORTRAN true CACHE BOOL "Enable generic FORTRAN code" )
  
  else()
  
    message( STATUS "HAXX Found a suitable SIMD instruction set -- Enabling Optimized code" )

    set( ENABLE_GENERIC_FORTRAN false CACHE BOOL "Enable generic FORTRAN code" )
  
  endif()


endif()
