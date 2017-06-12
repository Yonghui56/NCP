
#-------------------------------------------------------------------------------
# Scans the current directory and returns a list of subdirectories.
# Third parameter is 1 if you want relative paths returned.
# Usage: list_subdirectories(the_list_is_returned_here C:/cwd 1)
macro(list_subdirectories retval curdir return_relative)
#   MESSAGE (STATUS "curdir: ${curdir}")
#   file(GLOB sub-dir ${curdir}/*)
    file(GLOB sub-dir RELATIVE ${curdir} ${curdir}/*)
#   list(LENGTH sub-dir n_dir)
#   MESSAGE (STATUS "sub-dir (${n_dir}): ${sub-dir}")
    set(list_of_dirs "")
    foreach(dir ${sub-dir})
        if(IS_DIRECTORY ${curdir}/${dir})
            if (${return_relative})
                set(list_of_dirs ${list_of_dirs} ${dir})
            else()
                set(list_of_dirs ${list_of_dirs} ${curdir}/${dir})
            endif()
        endif()
    endforeach()
    set(${retval} ${list_of_dirs})
endmacro()

#-------------------------------------------------------------------------------

# look for directories in modules
#MESSAGE (STATUS "Look for subdirectory")
list_subdirectories(List_subdir ${CMAKE_CURRENT_SOURCE_DIR}/../UserModules 1)
LIST(LENGTH List_subdir n_dir)
MESSAGE (STATUS "${n_dir} processses are found.")

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../UserModules/InactiveModules.cmake)
LIST(LENGTH INACTIVE_MODULES n_inactive_dir)
MESSAGE (STATUS "${n_inactive_dir} processses are excluded. Inactive modules are defined in UserModules/CMakeLists.txt")
IF (n_inactive_dir GREATER 0)
    LIST(REMOVE_ITEM List_subdir ${INACTIVE_MODULES})
ENDIF()

FOREACH(dir ${List_subdir})
#    MESSAGE (STATUS "sub dir found: ${dir}")
    IF(EXISTS ../UserModules/${dir}/CMakeLists.txt)
        INCLUDE(../UserModules/${dir}/CMakeLists.txt)
    ENDIF()
    FILE(GLOB TempProcessLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ../UserModules/${dir}/*.h ../UserModules/${dir}/*.hpp ../UserModules/${dir}/*.tpp ../UserModules/${dir}/*.cpp ../UserModules/${dir}/*.txt ../UserModules/${dir}/*.cmake)
    SOURCE_GROUP( ogs6THMC\\UserModules\\${dir} FILES ${TempProcessLib_Files} )
    LIST(APPEND THMC_Files ${TempProcessLib_Files}) 
ENDFOREACH()

# ./
FILE(GLOB TempProcessLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h *.cpp *.txt *.cmake)
SOURCE_GROUP( ogs6THMC FILES ${TempProcessLib_Files} )
LIST(APPEND THMC_Files ${TempProcessLib_Files}) 

# ./modules
FILE(GLOB TempProcessLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ../UserModules/*.h ../UserModules/*.cpp ../UserModules/*.txt ../UserModules/*.cmake)
SOURCE_GROUP( ogs6THMC\\UserModules FILES ${TempProcessLib_Files} )
LIST(APPEND THMC_Files ${TempProcessLib_Files}) 

# ./MaterialLib
INCLUDE( MaterialLib/CMakeLists.txt )

# ./PhysicsLib
INCLUDE( PhysicsLib/CMakeLists.txt )

# ./ChemLib
INCLUDE( ChemLib/CMakeLists.txt )

# ./ProcessLib
INCLUDE( ProcessLib/CMakeLists.txt )

#
ADD_LIBRARY( ogs6THMC STATIC 
    ${MaterialLib_Files}
    ${ProcessLib_Files}
    ${ChemLib_Files}
    ${PhysicsLib_Files}
    ${THMC_Files}
)

target_link_libraries( ogs6THMC
    ogsCore
    ogsNumerics
    ogsFileIO
    tinyxml2
    muparser
)

SET_TARGET_PROPERTIES(ogs6THMC PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        .
        ..
        ../ogsCore
        ../ogsNumerics
        ../ogsFileIO
        ../external
        ../external/logog/include
        ../external/Eigen
        ../external/tinyxml2
		../external/muparser/include
)

IF(USE_EXTERNAL_LIS)
    INCLUDE_DIRECTORIES( ${LIS_INCLUDE_DIR} )
ELSE()
	IF (USE_LIS)
		INCLUDE_DIRECTORIES( ../external/lis/include )
	ENDIF()
ENDIF ()

IF(USE_PETSC)
        TARGET_LINK_LIBRARIES( ogs6THMC  ${PETSC_LIBRARIES} )
    IF (MSVC)
        SET_TARGET_PROPERTIES(ogs6THMC PROPERTIES
                LINK_FLAGS "
            /NODEFAULTLIB:\"libcmt.lib\"
                ")
    ENDIF(MSVC)
ENDIF(USE_PETSC)


