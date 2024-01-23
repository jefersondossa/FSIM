#ifndef PETSC_SOLVER_H
#define PETSC_SOLVER_H
#include "LinearSolver.h"

class PETScSolver : public LinearSolver
{
private:
    /* data */
public:
    PETScSolver(/* args */);
    ~PETScSolver();
};



#endif