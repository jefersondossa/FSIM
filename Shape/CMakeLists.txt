
file(GLOB shape 
    ShapeHexahedron.h
    ShapeHexahedron.cpp
    ShapeOneD.h
    ShapeOneD.cpp
    ShapePoint.h
    ShapePoint.cpp
    ShapeQuadrilateralLin.h
    ShapeQuadrilateralLin.cpp
    ShapeTetrahedronCub.h
    ShapeTetrahedronCub.cpp
    ShapeTetrahedronQua.h
    ShapeTetrahedronQua.cpp
    ShapeTetrahedronLin.h
    ShapeTetrahedronLin.cpp
    ShapeTriangleCub.h
    ShapeTriangleCub.cpp
    ShapeTriangleQua.h
    ShapeTriangleQua.cpp
    ShapeTriangleLin.h
    ShapeTriangleLin.cpp
    BoundaryShapeFunction.h
    BoundaryShapeFunction.cpp
)

set(SHAPE ${shape} PARENT_SCOPE)


