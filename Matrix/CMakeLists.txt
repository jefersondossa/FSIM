
file(GLOB matrix 
    DenseEigen.h
    PETScDataTypes.h
    SparseEigen.h
    MatrixType.h
    MatrixType.cpp
    PETScMatrix.h
    PETScMatrix.cpp
    EigenSpMatrix.h
    EigenSpMatrix.cpp
    Tensor.h
    Tensor.cpp
)

set(MATRIX ${matrix} PARENT_SCOPE)


