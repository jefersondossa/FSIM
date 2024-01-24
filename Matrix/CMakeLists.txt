
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
    LinearSolver.h
    LinearSolver.cpp
    PETScSolver.h
    PETScSolver.cpp
    EigenLinearSolver.h
    EigenLinearSolver.cpp
)

set(MATRIX ${matrix} PARENT_SCOPE)


