# Source files
FILE(GLOB BaseLib_Files BaseLib/*.h BaseLib/*.h.in BaseLib/*.cpp)

#SET( BaseLib_Files ${BaseLib_Files} PARENT_SCOPE )
SOURCE_GROUP( ogsCore\\BaseLib FILES ${BaseLib_Files} )

#MESSAGE (STATUS "BaseLib files are  ${BaseLib_Files}")

## Create the library
#ADD_LIBRARY( BaseLib STATIC ${Base_Files})

#SET_TARGET_PROPERTIES(BaseLib PROPERTIES LINKER_LANGUAGE CXX)

#INCLUDE_DIRECTORIES(
#        .
#)
