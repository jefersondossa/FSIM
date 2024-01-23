
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
    GlobalMatrix.h
    GlobalMatrix.cpp
    GlobalMatrixT.h
    GlobalMatrixT.cpp
    LinearSolver.h
    PETScSolver.h
    PETScSolver.cpp
    EigenLinearSolver.h
    EigenLinearSolver.cpp
)

set(MATRIX ${matrix} PARENT_SCOPE)


