#ifndef NONLINEARANALYSIS_H
#define NONLINEARANALYSIS_H

#include "LinearAnalysis.h"

class NonLinearAnalysis : public LinearAnalysis
{
protected:
    int fMaxIterations;
    double fTolerance;

public:
    NonLinearAnalysis() : LinearAnalysis(){};
    NonLinearAnalysis(CompMesh *cmesh, SolverType stype, double tol = 1.e-6, int maxIter = 10) : LinearAnalysis(cmesh,stype), fMaxIterations(maxIter), fTolerance(tol) {};
    NonLinearAnalysis(Arlequin* arl, SolverType stype, double tol = 1.e-6, int maxIter = 10) : LinearAnalysis(arl,stype), fMaxIterations(maxIter), fTolerance(tol) {};

    void UpdateSolution() override;

    void Run() override;
};

#endif