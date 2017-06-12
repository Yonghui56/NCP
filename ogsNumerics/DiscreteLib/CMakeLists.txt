#Source files grouped by a directory
# ./
FILE(GLOB DiscreteLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/*.h DiscreteLib/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib FILES ${DiscreteLib_Files} )
# ./Core
FILE(GLOB DiscreteLib_Core_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/Core/*.h DiscreteLib/Core/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\Core FILES ${DiscreteLib_Core_Files} )
# ./Utils
FILE(GLOB DiscreteLib_Utils_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/Utils/*.h DiscreteLib/Utils/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\Utils FILES ${DiscreteLib_Utils_Files} )
# ./DDC
FILE(GLOB DiscreteLib_DDC_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/DDC/*.h DiscreteLib/DDC/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\DDC FILES ${DiscreteLib_DDC_Files} )
# ./Serial
FILE(GLOB DiscreteLib_Searial_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/Serial/*.h DiscreteLib/Serial/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\Serial FILES ${DiscreteLib_Searial_Files} )
# ./OpenMP
FILE(GLOB DiscreteLib_OMP_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/OpenMP/*.h DiscreteLib/OpenMP/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\OpenMP FILES ${DiscreteLib_OMP_Files} )
# ./SerialNodeDdc
FILE(GLOB DiscreteLib_SearialNodeDdc_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/SerialNodeDdc/*.h DiscreteLib/SerialNodeDdc/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\SerialNodeDdc FILES ${DiscreteLib_SearialNodeDdc_Files} )
# ./lis
FILE(GLOB DiscreteLib_LIS_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/lis/*.h DiscreteLib/lis/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\lis FILES ${DiscreteLib_LIS_Files} )
IF(NOT USE_LIS OR NOT USE_MPI)
    LIST(REMOVE_ITEM DiscreteLib_LIS_Files 
        DiscreteLib/lis/LisMPIDiscreteVector.h 
        DiscreteLib/lis/LisMPIDiscreteVector.cpp
        DiscreteLib/lis/LisMPILinearSystem.h 
        DiscreteLib/lis/LisDiscreteSystem.h
        DiscreteLib/lis/LisDiscreteSystem.cpp
        )
ENDIF()

#
LIST(APPEND DiscreteLib_Files
    ${DiscreteLib_Core_Files} 
    ${DiscreteLib_Utils_Files} 
    ${DiscreteLib_DDC_Files} 
    ${DiscreteLib_Searial_Files} 
    ${DiscreteLib_OMP_Files} 
    ${DiscreteLib_SearialNodeDdc_Files} 
    ${DiscreteLib_LIS_Files} 
)


#ADD_LIBRARY( DiscreteLib STATIC 
#	${DiscreteLib_Files} 
#    ${DiscreteLib_Core_Files} 
#    ${DiscreteLib_LinearEquation_Files}
#    ${DiscreteLib_EquationId_Files} 
#    ${DiscreteLib_Assembler_Files} 
#    ${DiscreteLib_Utils_Files} 
#    ${DiscreteLib_DDC_Files} 
#    ${DiscreteLib_OpenMP_Files} 
#	${DiscreteLib_OGS5_Files} 
#)
#
#INCLUDE_DIRECTORIES(
#	.
#	..
#	../../ogsCore
#)
#
#IF (USE_MPI)
#	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#ENDIF (USE_MPI)
#
#target_link_libraries (
#	DiscreteLib
#	BaseLib
#        MathLib
#        MeshLib
#)
#
#SET_TARGET_PROPERTIES(DiscreteLib PROPERTIES LINKER_LANGUAGE CXX)
#