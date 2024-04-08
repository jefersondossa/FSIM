#ifndef PETSC_SOLVER_H
#define PETSC_SOLVER_H
#include "LinearSolver.h"

#ifdef HAS_PETSC

#include <petscksp.h> 

class PETScSolver : public LinearSolver
{
private:
    KSP               ksp;
    PC                pc;

public:
    PETScSolver(Analysis *an);
    virtual ~PETScSolver();

    void Solve() override;

    KSP *KSPSolver(){return &ksp;}
    PC  *Preconditioner(){return &pc;}
};

#endif

#endif