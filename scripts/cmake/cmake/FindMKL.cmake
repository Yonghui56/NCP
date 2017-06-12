# Find Intel Math Karnel Library

#if (MKL_INCLUDES AND MKL_LIBRARIES)
#  # Already in cache, be silent
#  set (MKL_FIND_QUIETLY TRUE)
#endif (MKL_INCLUDES AND MKL_LIBRARIES)

SET(MKL_DIR
    "${MKL_DIR}"
    CACHE
    PATH
    "Directory to search for MKL libraries")
IF (NOT MKL_DIR)
    MESSAGE (FATAL_ERROR "Please set MKL root diretory path" )
ENDIF()

FIND_PATH(MKL_INCLUDES NAMES mkl.h
    HINTS ${MKL_DIR} PATH_SUFFIXES include
)

# Tell if the unix system is on 64-bit base
if (${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
#if(CMAKE_SIZEOF_VOID_P MATCHES "8")
    MESSAGE (STATUS  "Look for 64-bit MKL" )
    find_library(MKL_LIB_CORE mkl_core 
      HINTS ${MKL_DIR} 
      PATH_SUFFIXES
      lib/intel64
      lib/em64t
    )

    find_library(MKL_LIB_INTEL mkl_intel_lp64 
      HINTS ${MKL_DIR} 
      PATH_SUFFIXES
      lib/intel64
      lib/em64t
    )

    find_library(MKL_LIB_SOLVER mkl_solver_lp64 
      HINTS ${MKL_DIR} 
      PATH_SUFFIXES
      lib/intel64
      lib/em64t
    )
    
    if (USE_OPENMP)
        if(CMAKE_C_COMPILER EQUAL "icc")
            find_library(MKL_LIB_INTEL_THREAD mkl_intel_thread 
              HINTS ${MKL_DIR} 
              PATH_SUFFIXES
              lib/intel64
              lib/em64t
            )
    
            find_library(MKL_LIB_IMOP5 iomp5 
              HINTS ${MKL_DIR} 
              PATH_SUFFIXES
              lib/intel64
              lib/em64t
            )
     
            find_library(MKL_LIB_PTHREAD pthread 
              HINTS ${MKL_DIR} 
              PATH_SUFFIXES
              lib/intel64
              lib/em64t
            )
            
            SET (MKL_LIB_THREAD 
                    "${MKL_LIB_INTEL_THREAD}" 
                    "${MKL_LIB_IMOP5}"
                    "${MKL_LIB_PTHREAD}")
        else()
            find_library(MKL_LIB_THREAD mkl_gnu_thread 
              HINTS ${MKL_DIR} 
              PATH_SUFFIXES
              lib/intel64
              lib/em64t
            )
        endif()
    else(USE_OPENMP)
        find_library(MKL_LIB_THREAD mkl_sequential
          HINTS ${MKL_DIR} 
          PATH_SUFFIXES
          lib/intel64
          lib/em64t
        )
    endif(USE_OPENMP)
    
    
else()
    MESSAGE (STATUS  "Look for 32-bit MKL" )
    
    find_library(MKL_LIB_CORE mkl_core 
      HINTS ${MKL_DIR} 
      PATH_SUFFIXES
      lib/32
    )

    find_library(MKL_LIB_INTEL mkl_intel 
      HINTS ${MKL_DIR} 
      PATH_SUFFIXES
      lib/32
    )

    find_library(MKL_LIB_SOLVER mkl_solver 
      HINTS ${MKL_DIR} 
      PATH_SUFFIXES
      lib/32
    )

    if (USE_OPENMP)
        if(CMAKE_C_COMPILER EQUAL "icc")
            find_library(MKL_LIB_INTEL_THREAD mkl_intel_thread 
              HINTS ${MKL_DIR} PATH_SUFFIXES lib/32
            )
    
            find_library(MKL_LIB_IMOP5 iomp5 
              HINTS ${MKL_DIR} PATH_SUFFIXES lib/32
            )
     
            find_library(MKL_LIB_PTHREAD pthread 
              HINTS ${MKL_DIR} PATH_SUFFIXES lib/32
            )
            
            SET (MKL_LIB_THREAD 
                    "${MKL_LIB_INTEL_THREAD}" 
                    "${MKL_LIB_IMOP5}"
                    "${MKL_LIB_PTHREAD}")
        else()
            find_library(MKL_LIB_THREAD mkl_gnu_thread 
              HINTS ${MKL_DIR} PATH_SUFFIXES lib/32
            )
        endif()
    else(USE_OPENMP)
        find_library(MKL_LIB_THREAD mkl_sequential
          HINTS ${MKL_DIR} PATH_SUFFIXES lib/32
        )
    endif(USE_OPENMP)

endif()

SET (MKL_LIBRARIES 
    "${MKL_LIB_INTEL}" 
#            "${MKL_LIB_SOLVER}"
    "${MKL_LIB_THREAD}"
    "${MKL_LIB_CORE}" 
        )
MESSAGE (STATUS "MKL libs: ${MKL_LIBRARIES}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDES MKL_LIBRARIES)
mark_as_advanced(MKL_INCLUDES MKL_LIBRARIES)

