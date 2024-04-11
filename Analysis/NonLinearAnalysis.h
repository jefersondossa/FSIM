#ifndef NONLINEARANALYSIS_H
#define NONLINEARANALYSIS_H

#include "LinearAnalysis.h"
#ifdef HAS_PETSC
#include <petscsnes.h>
#endif

enum NonLinearSolverType{ELinear,ENewtonRaphson,ELineSearch};

class NonLinearAnalysis : public LinearAnalysis
{
protected:
    int    fMaxIterations;
    double fTolerance;
    NonLinearSolverType fSolverType = ENewtonRaphson;
#ifdef HAS_PETSC
    SNES   fSNES;
#endif

public:
    NonLinearAnalysis() : LinearAnalysis(){};
    NonLinearAnalysis(CompMesh *cmesh, SolverType stype, double tol = 1.e-6, int maxIter = 10) : LinearAnalysis(cmesh,stype), fMaxIterations(maxIter), fTolerance(tol) {};
    NonLinearAnalysis(Arlequin* arl, SolverType stype, double tol = 1.e-6, int maxIter = 10) : LinearAnalysis(arl,stype), fMaxIterations(maxIter), fTolerance(tol) {};

    void SetMaxIter(int maxiter){
        fMaxIterations = maxiter;
    }

    void UpdateSolution() override;

    void Run() override;

#ifdef HAS_PETSC
    static PetscErrorCode FormFunction(SNES snes, Vec u,Vec b, void *ptr);
    static PetscErrorCode FormJacobian(SNES snes,Vec u,Mat A,Mat B,void *ptr);
#endif

};

#endif