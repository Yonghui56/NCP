#Source files grouped by a directory

FILE(GLOB ChemLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ChemLib/*.h ChemLib/*.cpp)
SOURCE_GROUP( ogsChemistry FILES ${ChemLib_Files} )

