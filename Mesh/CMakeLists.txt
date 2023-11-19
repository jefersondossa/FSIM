
file(GLOB mesh 
    Arlequin.h
    Arlequin.cpp
    CompMesh.h
    CompMesh.cpp
    Element.h
    Element.cpp
    ElementT.h
    ElementT.cpp
    FSInteraction.h
    FSInteraction.cpp
    Node.h
    Node.cpp
    IntPointData.h
    ElementWithMem.h
    ElementWithMem.cpp
)

set(MESH ${mesh} PARENT_SCOPE)


