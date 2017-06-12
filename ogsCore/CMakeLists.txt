
# include subdirectories
INCLUDE( BaseLib/CMakeLists.txt )
INCLUDE( MathLib/CMakeLists.txt )
INCLUDE( GeoLib/CMakeLists.txt )
INCLUDE( MeshLib/CMakeLists.txt )

#MESSAGE (STATUS "GeoLib files are ${GeoLib_Files}")
#MESSAGE (STATUS "MathLib files are ${MathLib_Files}")

ADD_LIBRARY( ogsCore STATIC 
	${BaseLib_Files}
	${MathLib_Files}
    ${GeoLib_Files}
    ${MeshLib_Files}
)

SET_TARGET_PROPERTIES(ogsCore PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        .
        ../external
        ../external/logog/include
        ../external/Eigen
        ../external/tinyxml2
)

TARGET_LINK_LIBRARIES( ogsCore 
    tinyxml2
    logog 
    )

IF(USE_EXTERNAL_LIS)
    INCLUDE_DIRECTORIES( ${LIS_INCLUDE_DIR} )
    TARGET_LINK_LIBRARIES( ogsCore ${LIS_LIBRARIES} )
ELSE()
	IF (USE_LIS)
		INCLUDE_DIRECTORIES( ../external/lis/include )
		TARGET_LINK_LIBRARIES( ogsCore lis )
	ENDIF()
ENDIF ()

IF(USE_PETSC)
    TARGET_LINK_LIBRARIES( ogsCore  ${PETSC_LIBRARIES} )
    IF (MSVC)
        SET_TARGET_PROPERTIES(ogsCore PROPERTIES
                LINK_FLAGS "
            /NODEFAULTLIB:\"libcmt.lib\"
                ")
    ENDIF( MSVC)
ENDIF(USE_PETSC)

IF(USE_MKL)
    INCLUDE_DIRECTORIES( ${MKL_INCLUDES} )
    TARGET_LINK_LIBRARIES( ogsCore ${MKL_LIBRARIES} )
ENDIF ()
