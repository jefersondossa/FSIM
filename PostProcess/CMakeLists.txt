
file(GLOB postprocess 
    VTUGenerator.h
    VTUGenerator.cpp
    XDMFGenerator.h
    XDMFGenerator.cpp
    GraphMesh.h
    GraphMesh.cpp
)

set(POSTPROCESS ${postprocess} PARENT_SCOPE)


