
file(GLOB preprocess 
    CompMeshTools.h
    CompMeshTools.cpp
    GmshTools.h
    GmshTools.cpp
    OpenFOAMWriter.h
    OpenFOAMWriter.cpp
    # ProblemParameters.h
)

set(PREPROCESS ${preprocess} PARENT_SCOPE)


