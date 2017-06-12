#Source files grouped by a directory
#.
FILE(GLOB FemLib_Files FemLib/*.h FemLib/*.cpp)
SOURCE_GROUP( ogsNumerics\\FemLib FILES ${FemLib_Files})
#Core
FILE(GLOB FemLib_Core_Files FemLib/Core/*.h FemLib/Core/*.cpp)
SOURCE_GROUP( ogsNumerics\\FemLib\\Core FILES ${FemLib_Core_Files})
#Core/Element
FILE(GLOB FemLib_Core_Element_Files FemLib/Core/Element/*.h FemLib/Core/Element/*.cpp)
SOURCE_GROUP( ogsNumerics\\FemLib\\Core\\Element FILES ${FemLib_Core_Element_Files})
#Core/CoordinatesMapping
FILE(GLOB FemLib_Core_CoordinatesMapping_Files FemLib/Core/CoordinatesMapping/*.h FemLib/Core/CoordinatesMapping/*.cpp)
SOURCE_GROUP( ogsNumerics\\FemLib\\Core\\CoordinatesMapping FILES ${FemLib_Core_CoordinatesMapping_Files})
#Core/Shape
FILE(GLOB FemLib_Core_Shape_Files FemLib/Core/Shape/*.h FemLib/Core/Shape/*.cpp)
SOURCE_GROUP( ogsNumerics\\FemLib\\Core\\Shape FILES ${FemLib_Core_Shape_Files})
#Core/Extrapolation
FILE(GLOB FemLib_Core_Extrapolation_Files FemLib/Core/Extrapolation/*.h FemLib/Core/Extrapolation/*.cpp)
SOURCE_GROUP( ogsNumerics\\FemLib\\Core\\Extrapolation FILES ${FemLib_Core_Extrapolation_Files})
#Core/Integration
FILE(GLOB FemLib_Core_Integration_Files FemLib/Core/Integration/*.h FemLib/Core/Integration/*.cpp)
SOURCE_GROUP( ogsNumerics\\FemLib\\Core\\Integration FILES ${FemLib_Core_Integration_Files})
#BC
FILE(GLOB FemLib_BC_Files FemLib/BC/*.h FemLib/BC/*.cpp)
SOURCE_GROUP( ogsNumerics\\FemLib\\BC FILES ${FemLib_BC_Files})
#Function
FILE(GLOB FemLib_Function_Files FemLib/Function/*.h FemLib/Function/*.cpp)
SOURCE_GROUP( ogsNumerics\\FemLib\\Function FILES ${FemLib_Function_Files})
#Post
FILE(GLOB FemLib_Post_Files FemLib/Post/*.h FemLib/Post/*.cpp)
SOURCE_GROUP( ogsNumerics\\FemLib\\Post FILES ${FemLib_Post_Files})
#Tools
FILE(GLOB FemLib_Tools_Files FemLib/Tools/*.h FemLib/Tools/*.cpp)
SOURCE_GROUP( ogsNumerics\\FemLib\\Tools FILES ${FemLib_Tools_Files})

#
LIST(APPEND FemLib_Files
	${FemLib_Files}
	${FemLib_Core_Files}
	${FemLib_Core_Element_Files}
    ${FemLib_Core_Shape_Files}
    ${FemLib_Core_Extrapolation_Files}
    ${FemLib_Core_Integration_Files}
    ${FemLib_Core_CoordinatesMapping_Files}
	${FemLib_BC_Files}
	${FemLib_Function_Files}
	${FemLib_Post_Files}
    ${FemLib_Tools_Files}
)

#ADD_LIBRARY(FemLib STATIC
#	${FemLib_Files}
#	${FemLib_Core_Files}
#	${FemLib_Core_Element_Files}
#    ${FemLib_Core_Shape_Files}
#    ${FemLib_Core_Extrapolation_Files}
#    ${FemLib_Core_Integration_Files}
#    ${FemLib_Core_CoordinatesMapping_Files}
#	${FemLib_BC_Files}
#	${FemLib_Function_Files}
#	${FemLib_Post_Files}
#    ${FemLib_Tools_Files}
#)
#
#INCLUDE_DIRECTORIES(
#	.
#	..
#	../../ogsCore
#)
#
#
#IF (PARALLEL_USE_MPI)
#	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#ENDIF (PARALLEL_USE_MPI)
#
#target_link_libraries (FemLib
#	BaseLib
#	MathLib
#	GeoLib
#	MeshLib
#	DiscreteLib
#)
#
#SET_TARGET_PROPERTIES(FemLib PROPERTIES LINKER_LANGUAGE CXX)
#