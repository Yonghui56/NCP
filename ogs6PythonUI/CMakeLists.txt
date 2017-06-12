#Source files grouped by a directory
FILE(GLOB ogsPython_Files 
	example.i 
	*.h *.cpp)
SOURCE_GROUP( ogsPython FILES ${ogsPython_Files} )

INCLUDE_DIRECTORIES(
	.
	..
)



#SET(CMAKE_SWIG_FLAGS "")
SET(CMAKE_SWIG_FLAGS -shadow)
SET(CMAKE_SWIG_OUTDIR ${LIBRARY_OUTPUT_PATH}/${CMAKE_CFG_INTDIR})

SET_SOURCE_FILES_PROPERTIES(example.i PROPERTIES CPLUSPLUS ON)
#SET_SOURCE_FILES_PROPERTIES(example.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_MODULE( ogsPython python 
	${ogsPython_Files}
	)
SWIG_LINK_LIBRARIES(ogsPython
	BaseLib
	MathLib
	GeoLib
	MeshLib
	DiscreteLib
	FemLib
	NumLib
	 ${PYTHON_LIBRARIES}
	 )

IF (LIS_FOUND)
	SWIG_LINK_LIBRARIES( ogsPython ${LIS_LIBRARIES} )
ENDIF (LIS_FOUND)

#SET_TARGET_PROPERTIES(ogsPython PROPERTIES LINKER_LANGUAGE CXX)
