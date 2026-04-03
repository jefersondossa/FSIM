
file(GLOB postprocess 
    VTUGenerator.h
    VTUGenerator.cpp
    XDMFGenerator.h
    XDMFGenerator.cpp
    GraphMesh.h
    GraphMesh.cpp
)

set(POSTPROCESS ${postprocess} PARENT_SCOPE)

add_library(vtu_reader VTUReader.cpp)
target_include_directories(vtu_reader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

