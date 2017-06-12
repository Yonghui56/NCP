# Source files
FILE(GLOB FEMIO_FILES FemIO/*.h FemIO/*.cpp FemIO/ogs5/*.h FemIO/ogs5/*.cpp)
SOURCE_GROUP(ogsFileIO\\FemIO\\ogs5 FILES ${FEMIO_FILES})


#
#LIST(APPEND MeshLib_Files
#	${MeshLib_Core_Files} 
#)

