cmake_minimum_required (VERSION 2.8.4)
enable_testing()
project (muparser)
include(CTest)
# muparser needs thread support on linux
# find_package( Threads )
set( RUNTIME_OUTPUT_DIRECTORY bin/ )
set( ARCHIVE_OUTPUT_DIRECTORY bin/ )
set( LIBRARY_OUTPUT_DIRECTORY bin/ )
set( CMAKE_BINARY_DIR build/ )
set( CMAKE_LEGACY_CYGWIN_WIN32 0 )
if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-variadic-macros -Wall -Wno-long-long")
endif()
include_directories( include )
add_library( muparser
	src/muParser.cpp 
	src/muParserBase.cpp 
	src/muParserBytecode.cpp
    src/muParserCallback.cpp
    src/muParserDLL.cpp
    src/muParserError.cpp
    src/muParserInt.cpp
    src/muParserTest.cpp
    src/muParserTokenReader.cpp 	
)
set_target_properties(muparser PROPERTIES DEBUG_POSTFIX "d")
# add_executable( test-muparser src/muParserTest.cpp )
target_link_libraries( muparser ${CMAKE_THREAD_LIBS_INIT})
# add_test( NAME test-harness COMMAND test-muparser )
install(TARGETS muparser ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(DIRECTORY include/ DESTINATION "${CMAKE_INSTALL_PREFIX}/include/muparser"
          FILES_MATCHING PATTERN "*.h")

# Docs generation with Doxygen
# find_package( Doxygen )
# if( DOXYGEN_FOUND )
# 	add_custom_target (muparser-doc ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/doxyfile
# 		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#		COMMENT "Generating source code documentation with Doxygen." VERBATIM)
# endif() # DOXYGEN_FOUND