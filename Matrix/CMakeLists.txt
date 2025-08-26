
file(GLOB matrix 
    DenseEigen.h
    DenseEigen.cpp
    PETScDataTypes.h
    SparseEigen.h
    MatrixType.h
    MatrixType.cpp
    PETScMatrix.h
    PETScMatrix.cpp
    EigenSpMatrix.h
    EigenSpMatrix.cpp
    Tensor3D.h
    Tensor3D.cpp
    ArlequinMatRedEigen.h
    ArlequinMatRedEigen.cpp
    ArlequinMatRedPETSc.h
    ArlequinMatRedPETSc.cpp
)

set(MATRIX ${matrix} PARENT_SCOPE)


