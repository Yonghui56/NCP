#Source files grouped by a directory

#./
FILE(GLOB NumLib_Files 
	NumLib/*.h NumLib/*.cpp)
SOURCE_GROUP( ogsNumerics\\NumLib FILES ${NumLib_Files} )
#TimeStepping/
FILE(GLOB NumLib_TimeStepping_Files 
	NumLib/TimeStepping/*.h NumLib/TimeStepping/*.cpp)
SOURCE_GROUP( ogsNumerics\\NumLib\\TimeStepping FILES ${NumLib_TimeStepping_Files} )
#TransientCoupling/
FILE(GLOB NumLib_TransientCoupling_Files 
	NumLib/TransientCoupling/*.h NumLib/TransientCoupling/*.cpp)
SOURCE_GROUP( ogsNumerics\\NumLib\\TransientCoupling FILES ${NumLib_TransientCoupling_Files} )
#TransientAssembler/
FILE(GLOB NumLib_TransientAssembler_Files 
    NumLib/TransientAssembler/*.h NumLib/TransientAssembler/*.cpp)
SOURCE_GROUP( ogsNumerics\\NumLib\\TransientAssembler FILES ${NumLib_TransientAssembler_Files} )
#Nonlinear/
FILE(GLOB NumLib_Nonlinear_Files 
    NumLib/Nonlinear/*.h NumLib/Nonlinear/*.cpp)
SOURCE_GROUP( ogsNumerics\\NumLib\\Nonlinear FILES ${NumLib_Nonlinear_Files} )
#Function/
FILE(GLOB NumLib_Function_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	NumLib/Function/*.h 
	NumLib/Function/*.cpp)
SOURCE_GROUP( ogsNumerics\\NumLib\\Function FILES ${NumLib_Function_Files})
#IOSystem/
FILE(GLOB NumLib_IOSystem_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    NumLib/IOSystem/*.h 
    NumLib/IOSystem/*.cpp)
SOURCE_GROUP( ogsNumerics\\NumLib\\IOSystem FILES ${NumLib_IOSystem_Files})
#Coupling/
FILE(GLOB NumLib_Coupling_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    NumLib/Coupling/*.h 
    NumLib/Coupling/*.cpp)
SOURCE_GROUP( ogsNumerics\\NumLib\\Coupling FILES ${NumLib_Coupling_Files})
#Coupling/Algorithm
FILE(GLOB NumLib_Coupling_Algorithm_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    NumLib/Coupling/Algorithm/*.h 
    NumLib/Coupling/Algorithm/*.cpp)
SOURCE_GROUP( ogsNumerics\\NumLib\\Coupling\\Algorithm FILES ${NumLib_Coupling_Algorithm_Files})
#SystemOfEquations/
FILE(GLOB NumLib_SystemOfEquations_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    NumLib/SystemOfEquations/*.h 
    NumLib/SystemOfEquations/*.cpp)
SOURCE_GROUP( ogsNumerics\\NumLib\\SystemOfEquations FILES ${NumLib_SystemOfEquations_Files})

#
LIST(APPEND NumLib_Files
    ${NumLib_Nonlinear_Files} 
	${NumLib_TimeStepping_Files} 
    ${NumLib_TransientCoupling_Files} 
	${NumLib_TransientAssembler_Files} 
    ${NumLib_Coupling_Files}
    ${NumLib_Coupling_Algorithm_Files}
    ${NumLib_Function_Files}
    ${NumLib_IOSystem_Files}
    ${NumLib_SystemOfEquations_Files}
)

#ADD_LIBRARY( NumLib STATIC 
#	${NumLib_Files} 
##	${NumLib_Core_Files} 
#    ${NumLib_Nonlinear_Files} 
#    ${NumLib_FileIO_Files} 
#	${NumLib_TimeStepping_Files} 
#    ${NumLib_TransientCoupling_Files} 
#	${NumLib_TransientAssembler_Files} 
#    ${NumLib_Coupling_Files}
#    ${NumLib_Coupling_Algorithm_Files}
#    ${NumLib_Function_Files}
#    ${NumLib_IOSystem_Files}
#    ${NumLib_SystemOfEquations_Files}
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
#	NumLib
#	BaseLib
#        NumLib
#        MeshLib
#	DiscreteLib
#)
#SET_TARGET_PROPERTIES(NumLib PROPERTIES LINKER_LANGUAGE CXX)
#