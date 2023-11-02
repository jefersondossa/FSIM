
file(GLOB mesh 
    Arlequin.h
    Arlequin.cpp
    Boundary.h
    CompMesh.h
    ElBoundary.cpp
    Element.h
    Element.cpp
    ElementT.h
    ElementT.cpp
    FSInteraction.h
    FSInteraction.cpp
    Glue.h
    Node.h
    Node.cpp
    IntPointData.h
    Fluid.h 
    Fluid.cpp
)

set(MESH ${mesh} PARENT_SCOPE)


