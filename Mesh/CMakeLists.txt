
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
    Connect.h
    Connect.cpp 
    MixedCompMesh.h
    MixedCompMesh.cpp
    GeoMesh.h
    GeoMesh.cpp
    GeoElement.h
    GeoElement.cpp
    GeoElementT.h
    GeoElementT.cpp
    ElementMixed.h
    ElementMixed.cpp
    ElementEnriched.h
    ElementEnriched.cpp
)

set(MESH ${mesh} PARENT_SCOPE)


