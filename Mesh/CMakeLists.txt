
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
    ElementTransient.h
    ElementTransient.cpp
    BoundaryCoupling.h
    BoundaryCoupling.cpp
    BndCondCoupling.h
    BndCondCoupling.cpp
    LagrangeMultCoupling.h
    LagrangeMultCoupling.cpp
    ElementLagrangeMultiplier.h
    ElementLagrangeMultiplier.cpp
)

set(MESH ${mesh} PARENT_SCOPE)


