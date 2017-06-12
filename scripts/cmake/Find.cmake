######################
### Find tools     ###
######################

## Unix tools ##
# Date
FIND_PROGRAM(DATE_TOOL_PATH date PATHS ${MSYSGIT_BIN_DIR})
# Grep
FIND_PROGRAM(GREP_TOOL_PATH grep PATHS ${MSYSGIT_BIN_DIR})
# Unzip
FIND_PROGRAM(UNZIP_TOOL_PATH unzip PATHS ${MSYSGIT_BIN_DIR})

# Find dot tool from graphviz
FIND_PROGRAM(DOT_TOOL_PATH dot DOC "Dot tool from graphviz")

# Find doxygen
FIND_PACKAGE(Doxygen)

# Find gnu profiler gprof
FIND_PROGRAM(GPROF_PATH gprof DOC "GNU profiler gprof")

FIND_PACKAGE(cppcheck)

# Find Git
FIND_PACKAGE(Git)

# msysGit on Windows
IF(WIN32 AND GIT_FOUND)
    FIND_PACKAGE(MsysGit)
ENDIF() # WIN32 AND GIT_FOUND

######################
### Find libraries ###
######################

# Clang does not have OpenMP support atm, see https://github.com/ufz/ogs/issues/8
#IF(NOT COMPILER_IS_CLANG)
#	FIND_PACKAGE(OpenMP)
#ENDIF () # !clang
#IF(OPENMP_FOUND)
#	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#ENDIF()

