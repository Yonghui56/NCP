INCLUDE(ResetConfigurations)        # To Debug, Release, RelWithDbgInfo
INCLUDE(SetDefaultBuildType)
INCLUDE(DisableCompilerFlag)
SET_DEFAULT_BUILD_TYPE(Debug)
INCLUDE(MSVCMultipleProcessCompile) # /MP switch (multi processor) for VS

# Set compiler helper variables
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(COMPILER_IS_CLANG TRUE)
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(COMPILER_IS_GCC TRUE)
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    SET(COMPILER_IS_INTEL TRUE)
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    SET(COMPILER_IS_MSVC TRUE)
ENDIF()

### GNU C/CXX compiler
IF(COMPILER_IS_GCC)
		get_gcc_version(GCC_VERSION)
        IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
                MESSAGE(STATUS "Set GCC release flags")
				IF(APPLE AND GCC_VERSION VERSION_LESS "4.3" AND NOT "${CMAKE_GENERATOR}" STREQUAL "Xcode" )
					# -march=native does not work here when on normal gcc compiler
					# see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=33144
					SET(CMAKE_CXX_FLAGS "-O3 -mtune=native -msse4.2 -DNDEBUG")
				ELSE()
                	SET(CMAKE_CXX_FLAGS "-O3 -march=native -mtune=native -msse4.2 -DNDEBUG")
				ENDIF()
        ENDIF()
        # -g
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wall -Wextra -fno-nonansi-builtins -std=gnu++0x")
        ADD_DEFINITIONS( -DGCC )
ENDIF() # COMPILER_IS_GCC

### Intel compiler
IF (COMPILER_IS_INTEL)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
                MESSAGE(STATUS "Set Intel release flags")
                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
        ENDIF()
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST -no-prec-div")
ENDIF() # COMPILER_IS_INTEL

### Clang compiler
IF (COMPILER_IS_CLANG)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
                MESSAGE(STATUS "Set Clang release flags")
                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
        ENDIF()
         # need -std=gnu89 to avoid linking errors of multiple definitions
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu89")
ENDIF() # COMPILER_IS_CLANG

# Profiling
IF (OGS_PROFILE)
	IF( NOT CMAKE_BUILD_TYPE STREQUAL "Release" )
		MESSAGE(STATUS "When using profiling you should set CMAKE_BUILD_TYPE to Release.")
	ENDIF()
	SET(PROFILE_FLAGS "-pg -fno-omit-frame-pointer -O2 -DNDEBUG")
	# clang compiler does not know the following flags
	IF(NOT COMPILER_IS_CLANG)
		SET(PROFILE_FLAGS "${PROFILE_FLAGS} -fno-inline-functions-called-once -fno-optimize-sibling-calls")
	ENDIF()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROFILE_FLAGS}")
ENDIF (OGS_PROFILE)

### Windows
IF (WIN32)
	## For Visual Studio compiler
	IF (MSVC)
		ADD_DEFINITIONS(
			-D_CRT_SECURE_NO_WARNINGS
			-D_CRT_NONSTDC_NO_WARNINGS
			-D_CRT_XNONSTDC_NO_WARNINGS
			-D__restrict__=__restrict   # this fixes #5
			-DNOMINMAX # This fixes compile errors with std::numeric_limits<T>::min() / max()
		)
		# Sets warning level 3 and ignores some warnings
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /wd4290 /wd4267")

		DisableCompilerFlag(DEBUG /RTC1)
	# cygwin
	ELSE (MSVC)
		MESSAGE (STATUS "Might be GCC under cygwin.")
		ADD_DEFINITIONS( -DGCC )
	ENDIF (MSVC)
ENDIF (WIN32)

# Missing OpenMP 3.0 implementation fix for Windows, this fixes #6
IF(MSVC)
	ADD_DEFINITIONS(-DOPENMP_LOOP_TYPE=int)
ELSE()
	ADD_DEFINITIONS(-DOPENMP_LOOP_TYPE=unsigned)
ENDIF()

# Use override and final keywords if a compiler supports them
SET(ENABLE_C11_OVERRIDE_FINAL FALSE)
IF(COMPILER_IS_GCC)
        get_gcc_version(GCC_VERSION)
        IF(GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
            SET(ENABLE_C11_OVERRIDE_FINAL TRUE)
        ENDIF()
ENDIF()
IF(COMPILER_IS_MSVC)
    SET(ENABLE_C11_OVERRIDE_FINAL TRUE)
ENDIF()
IF(ENABLE_C11_OVERRIDE_FINAL)
    MESSAGE(STATUS "ENABLE_C11_OVERRIDE_FINAL true")
    ADD_DEFINITIONS(-DOGS_ENABLE_DECL_OVERRIDE_FINAL)
ENDIF()
