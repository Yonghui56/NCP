#Source files grouped by a directory
# ./
FILE(GLOB PhysicsLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} PhysicsLib/*.h PhysicsLib/*.cpp)
SOURCE_GROUP( ogs6THMC\\PhysicsLib FILES ${PhysicsLib_Files} )
