# Add subdirectories with the projects
#ADD_SUBDIRECTORY( DiscreteLib )
#ADD_SUBDIRECTORY( FemLib )
#ADD_SUBDIRECTORY( MeshLib )
#ADD_SUBDIRECTORY( NumLib )
#ADD_SUBDIRECTORY( SolutionLib )
#INCLUDE( MeshLib/CMakeLists.txt )
INCLUDE( DiscreteLib/CMakeLists.txt )
INCLUDE( FemLib/CMakeLists.txt )
INCLUDE( NumLib/CMakeLists.txt )
INCLUDE( SolutionLib/CMakeLists.txt )

ADD_LIBRARY( ogsNumerics STATIC 
#	${MeshLib_Files}
	${DiscreteLib_Files}
	${FemLib_Files}
	${NumLib_Files}
	${SolutionLib_Files}
)

SET_TARGET_PROPERTIES(ogsNumerics PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        .
        ../ogsCore
        ../external
        ../external/logog/include
        ../external/Eigen
)

target_link_libraries( ogsNumerics
    ogsCore
)
