
# Source files
FILE(GLOB tinyxml2_Files *.h *.cpp)

#SET( BaseLib_Files ${BaseLib_Files} PARENT_SCOPE )
SOURCE_GROUP( tinyxml2 FILES ${tinyxml2_Files} )

INCLUDE_DIRECTORIES(
        .
)

ADD_LIBRARY( tinyxml2 STATIC 
	${tinyxml2_Files}
)

SET_TARGET_PROPERTIES(tinyxml2 PROPERTIES LINKER_LANGUAGE CXX)

