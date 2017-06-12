#Source files grouped by a directory
FILE(GLOB SolutionLib_Files SolutionLib/*.h SolutionLib/*.cpp)
SOURCE_GROUP( ogsNumerics\\SolutionLib FILES ${SolutionLib_Files} )
# ./Core
FILE(GLOB SolutionLib_Core_Files SolutionLib/Core/*.h SolutionLib/Core/*.cpp)
SOURCE_GROUP( ogsNumerics\\SolutionLib\\Core FILES ${SolutionLib_Core_Files} )
# ./Fem
FILE(GLOB SolutionLib_Fem_Files SolutionLib/Fem/*.h SolutionLib/Fem/*.cpp)
SOURCE_GROUP( ogsNumerics\\SolutionLib\\Fem FILES ${SolutionLib_Fem_Files} )

#
LIST(APPEND SolutionLib_Files
	${SolutionLib_Files}
    ${SolutionLib_Core_Files}
    ${SolutionLib_Fem_Files}
)
