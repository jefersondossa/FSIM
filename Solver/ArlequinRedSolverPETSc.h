#ifndef ArlequinRedSolverPETSc_H
#define ArlequinRedSolverPETSc_H
#include "LinearSolver.h"

#ifdef HAS_PETSC

#include <petscksp.h> 
#include "PETScSolver.h"
#include "ArlequinMatRedPETSc.h"

class ArlequinRedSolverPETSc : public LinearSolver
{
private:
    KSP               kspU0, kspU1, kspLambda;
    PC                pcU0, pcU1, pcLambda;

public:
    ArlequinRedSolverPETSc(Analysis *an);
    virtual ~ArlequinRedSolverPETSc();

    void Solve() override;

    // KSP *KSPSolver(){return &ksp;}
    // PC  *Preconditioner(){return &pc;}
};

#endif

#endif