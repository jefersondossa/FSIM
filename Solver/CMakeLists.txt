
file(GLOB solver 
    LinearSolver.h
    LinearSolver.cpp
    PETScSolver.h
    PETScSolver.cpp
    EigenLinearSolver.h
    EigenLinearSolver.cpp
    ArlequinRedSolverEigen.h
    ArlequinRedSolverEigen.cpp
    ArlequinRedSolverPETSc.h
    ArlequinRedSolverPETSc.cpp
)

set(SOLVER ${solver} PARENT_SCOPE)


