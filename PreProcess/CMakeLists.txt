
file(GLOB preprocess 
    CompMeshTools.h
    CompMeshTools.cpp
    GmshTools.h
    GmshTools.cpp
    OpenFOAMWriter.h
    OpenFOAMWriter.cpp
    GeoMeshTools.h
    GeoMeshTools.cpp
)

set(PREPROCESS ${preprocess} PARENT_SCOPE)


