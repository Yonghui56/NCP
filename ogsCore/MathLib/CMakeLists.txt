#Source files grouped by a directory
#./
FILE(GLOB MathLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	MathLib/*.h 
	MathLib/*.cpp)
SOURCE_GROUP( ogsCore\\MathLib FILES ${MathLib_Files} )
#LinAlg/
FILE(GLOB MathLib_LinAlg_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	MathLib/LinAlg/*.h 
	MathLib/LinAlg/*.cpp)
SOURCE_GROUP( ogsCore\\MathLib\\LinAlg FILES ${MathLib_LinAlg_Files} )
#LinAlg/Dense/
FILE(GLOB MathLib_LinAlg_Dense_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	MathLib/LinAlg/Dense/*.h 
	MathLib/LinAlg/Dense/*.cpp)
SOURCE_GROUP( ogsCore\\MathLib\\LinAlg\\Dense FILES ${MathLib_LinAlg_Dense_Files} )
#LinAlg/Sparse/
FILE(GLOB MathLib_LinAlg_Sparse_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	MathLib/LinAlg/Sparse/*.h 
	MathLib/LinAlg/Sparse/*.cpp)
IF (NOT USE_EIGEN)
	LIST(REMOVE_ITEM MathLib_LinAlg_Sparse_Files 
		MathLib/LinAlg/Sparse/EigenInterface.h)
ENDIF ()
SOURCE_GROUP( ogsCore\\MathLib\\LinAlg\\Sparse FILES ${MathLib_LinAlg_Sparse_Files} )
#LinAlg/Solvers/
FILE(GLOB MathLib_LinAlg_Solvers_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	MathLib/LinAlg/Solvers/*.h 
	MathLib/LinAlg/Solvers/*.cpp)
IF (NOT USE_BLAS_LAPACK)
	LIST(REMOVE_ITEM MathLib_LinAlg_Solvers_Files 
	        MathLib/LinAlg/Solvers/BiCGStab.h
	        MathLib/LinAlg/Solvers/BiCGStab.cpp
	        MathLib/LinAlg/Solvers/CG.h
	        MathLib/LinAlg/Solvers/CG.cpp
	)
ENDIF ()
SOURCE_GROUP( ogsCore\\MathLib\\LinAlg\\Solvers FILES ${MathLib_LinAlg_Solvers_Files} )
#LinAlg/Preconditioner/
FILE(GLOB MathLib_LinAlg_Preconditioner_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	MathLib/LinAlg/Preconditioner/*.h 
	MathLib/LinAlg/Preconditioner/*.cpp)
SOURCE_GROUP( ogsCore\\MathLib\\LinAlg\\Preconditioner FILES ${MathLib_LinAlg_Preconditioner_Files} )
#LinAlg/LinearEquation/
FILE(GLOB MathLib_LinAlg_LinearEquations_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	MathLib/LinAlg/LinearEquation/*.h 
	MathLib/LinAlg/LinearEquation/*.cpp)
IF (NOT USE_BLAS_LAPACK)
	LIST(REMOVE_ITEM MathLib_LinAlg_LinearEquations_Files 
		MathLib/LinAlg/LinearEquation/SparseLinearEquation.h
		MathLib/LinAlg/LinearEquation/SparseLinearEquation.cpp
	)
ENDIF ()
IF (NOT USE_LIS)
	LIST(REMOVE_ITEM MathLib_LinAlg_LinearEquations_Files 
		MathLib/LinAlg/LinearEquation/LisLinearEquation.h
		MathLib/LinAlg/LinearEquation/LisLinearEquation.cpp
	)
ENDIF ()
IF(NOT (USE_MPI AND USE_LIS))
	LIST(REMOVE_ITEM MathLib_LinAlg_LinearEquations_Files 
		MathLib/LinAlg/LinearEquation/LisMPILinearEquation.h
		MathLib/LinAlg/LinearEquation/LisMPILinearEquation.cpp
	)
ENDIF()

IF (NOT USE_PETSC)
	LIST(REMOVE_ITEM MathLib_LinAlg_LinearEquations_Files 
		MathLib/LinAlg/LinearEquation/PETScLinearEquation.h
		MathLib/LinAlg/LinearEquation/PETScLinearEquation.cpp
	)
ENDIF ()
IF(NOT (USE_MPI AND USE_PETSC))
    LIST(REMOVE_ITEM MathLib_LinAlg_LinearEquations_Files 
        MathLib/LinAlg/LinearEquation/PETScMPILinearEquation.h
        MathLib/LinAlg/LinearEquation/PETScMPILinearEquation.cpp
    )
ENDIF()

IF (NOT USE_PARDISO)
    MESSAGE (STATUS "USE_PARDISO is OFF. Exclude Paridso files from a project.")
    LIST(REMOVE_ITEM MathLib_LinAlg_LinearEquations_Files 
        MathLib/LinAlg/LinearEquation/PARDISOLinearEquation.h
        MathLib/LinAlg/LinearEquation/PARDISOLinearEquation.cpp
    )
ENDIF ()

SOURCE_GROUP( ogsCore\\MathLib\\LinAlg\\LinearEquation FILES ${MathLib_LinAlg_LinearEquations_Files})
#Nonlinear/
FILE(GLOB MathLib_Nonlinear_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	MathLib/Nonlinear/*.h 
	MathLib/Nonlinear/*.cpp)
SOURCE_GROUP( ogsCore\\MathLib\\Nonlinear FILES ${MathLib_Nonlinear_Files})
#Integration/
FILE(GLOB MathLib_Integration_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	MathLib/Integration/*.h 
	MathLib/Integration/*.cpp)
SOURCE_GROUP( ogsCore\\MathLib\\Integration FILES ${MathLib_Integration_Files})
#ODE/
FILE(GLOB MathLib_ODE_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	MathLib/ODE/*.h 
	MathLib/ODE/*.cpp)
SOURCE_GROUP( ogsCore\\MathLib\\ODE FILES ${MathLib_ODE_Files})
#Interpolation/
FILE(GLOB MathLib_Interpolation_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    MathLib/Interpolation/*.h 
    MathLib/Interpolation/*.cpp)
SOURCE_GROUP( ogsCore\\MathLib\\Interpolation FILES ${MathLib_Interpolation_Files})

#
LIST(APPEND MathLib_Files
    ${MathLib_Integration_Files}
    ${MathLib_Interpolation_Files}
	${MathLib_LinAlg_Files} 
	${MathLib_LinAlg_Dense_Files}
	${MathLib_LinAlg_Sparse_Files}
	${MathLib_LinAlg_Solvers_Files}
	${MathLib_LinAlg_Preconditioner_Files}
    ${MathLib_LinAlg_LinearEquations_Files}
	${MathLib_Nonlinear_Files}
	${MathLib_ODE_Files}
)

#FOREACH ( f ${MathLib_Files} )
#	LIST( APPEND MathLib_ABS_Files 
#		MathLib/${f})
#ENDFOREACH ( f )

#SET( MathLib_Files ${MathLib_ABS_Files} PARENT_SCOPE )

## Create the library
#ADD_LIBRARY( MathLib STATIC 
#	${MathLib_Files} 
#    ${MathLib_Integration_Files}
#    ${MathLib_Interpolation_Files}
#	${MathLib_LinAlg_Files} 
#	${MathLib_LinAlg_Dense_Files}
#	${MathLib_LinAlg_Sparse_Files}
#	${MathLib_LinAlg_Solvers_Files}
#	${MathLib_LinAlg_Preconditioner_Files}
#    ${MathLib_LinAlg_LinearEquations_Files}
#	${MathLib_Nonlinear_Files}
#	${MathLib_ODE_Files}
#)
#
##Include
#INCLUDE_DIRECTORIES (
#	.
#	..
#)
IF (USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (USE_MPI)

IF(USE_EXTERNAL_LIS)
	INCLUDE_DIRECTORIES( ${LIS_INCLUDE_DIR} )
ENDIF(USE_EXTERNAL_LIS)

#SET_TARGET_PROPERTIES(MathLib PROPERTIES LINKER_LANGUAGE CXX)
#
#