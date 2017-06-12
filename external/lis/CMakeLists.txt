#####################
### LIS ###
#####################

# Configurations
#OPTION(ENABLE_OPENMP "Use OpenMP" OFF)
#OPTION(ENABLE_MPI "Use MPI" OFF)
#OPTION(ENABLE_FORTRAN "build Fortran 77 interfaces" OFF)
#OPTION(ENABLE_SAAMG "build SA-AMG preconditioner" OFF)
#OPTION(ENABLE_QUAD "enable quadruple precision operations" OFF)
#OPTION(ENABLE_SSE2 "use Intel Streaming SIMD Extensions" ON)
#OPTION(BUILD_TEST "Build tests" OFF)

SET(ENABLE_OPENMP ${USE_OPENMP})
SET(ENABLE_MPI  ${USE_MPI})

###########################################################################
IF (MSVC)
	ADD_DEFINITIONS(-DHAVE_CONFIG_WIN32_H)
ELSE()
	ADD_DEFINITIONS(-DHAVE_CONFIG_H)
ENDIF()

######################
### Subdirectories ###
######################
# Add subdirectories with the projects
ADD_SUBDIRECTORY( src )
#IF (BUILD_TEST)
#	ADD_SUBDIRECTORY( test )
#ENDIF()
