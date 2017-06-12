# Source files
FILE(GLOB GeoLib_Files GeoLib/*.h GeoLib/*.cpp)
SOURCE_GROUP( ogsCore\\GeoLib FILES ${GeoLib_Files})


#SET( GeoLib_Files ${GeoLib_Files} PARENT_SCOPE )
#MESSAGE (STATUS "GeoLib files are ${GeoLib_Files}")


# Create the library
#ADD_LIBRARY(ogsGeometry STATIC 
#	${GeoLib_Files} 
#	${GeoLib_Core_Files}
#	${GeoLib_Shape_Files}
#)

#include_directories(
#	.
#	../ogsCore
#)


#target_link_libraries (
#	GeoLib
#	BaseLib
#	MathLib
#)

#SET_TARGET_PROPERTIES(ogsGeometry PROPERTIES LINKER_LANGUAGE CXX)

