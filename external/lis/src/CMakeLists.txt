
FILE(GLOB temp_Files esolver/*.h esolver/*.c)
SOURCE_GROUP( lis\\src\\esolver FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
IF (ENABLE_SAAMG OR ENABLE_FORTRAN)
	FILE(GLOB temp_Files fortran/*.h fortran/*.c fortran/*.f fortran/*.F)
	SOURCE_GROUP( lis\\src\\fortran FILES ${temp_Files} )
	LIST(APPEND lis_src_Files ${temp_Files}) 
ENDIF()
IF (ENABLE_SAAMG)
	FILE(GLOB temp_Files fortran/amg/*.h fortran/amg/*.c fortran/amg/*.f fortran/amg/*.F)
	SOURCE_GROUP( lis\\src\\fortran\\amg FILES ${temp_Files} )
	LIST(APPEND lis_src_Files ${temp_Files}) 
ENDIF()
FILE(GLOB temp_Files matrix/*.h matrix/*.c)
SOURCE_GROUP( lis\\src\\matrix FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
FILE(GLOB temp_Files matvec/*.h matvec/*.c)
SOURCE_GROUP( lis\\src\\matvec FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
FILE(GLOB temp_Files precision/*.h precision/*.c)
SOURCE_GROUP( lis\\src\\precision FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
FILE(GLOB temp_Files precon/*.h precon/*.c)
SOURCE_GROUP( lis\\src\\precon FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
FILE(GLOB temp_Files solver/*.h solver/*.c)
SOURCE_GROUP( lis\\src\\solver FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
FILE(GLOB temp_Files system/*.h system/*.c)
SOURCE_GROUP( lis\\src\\system FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
FILE(GLOB temp_Files vector/*.h vector/*.c)
SOURCE_GROUP( lis\\src\\vector FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 

FILE(GLOB lis_include_Files ../include/*.h)
SOURCE_GROUP( lis\\include FILES ${lis_include_Files} )

INCLUDE_DIRECTORIES(
	.
	../include
)

# Create the executable
ADD_LIBRARY( lis
        ${lis_include_Files}
        ${lis_src_Files}
)




IF (USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
        TARGET_LINK_LIBRARIES( tests ${MPI_C_LIBRARIES} )
ENDIF (USE_MPI)


SET_TARGET_PROPERTIES(lis PROPERTIES LINKER_LANGUAGE C)

