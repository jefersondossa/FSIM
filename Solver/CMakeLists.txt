
file(GLOB solver 
    LinearSolver.h
    LinearSolver.cpp
    PETScSolver.h
    PETScSolver.cpp
    EigenLinearSolver.h
    EigenLinearSolver.cpp
)

set(SOLVER ${solver} PARENT_SCOPE)


