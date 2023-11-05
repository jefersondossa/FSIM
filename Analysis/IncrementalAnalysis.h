#ifndef INCREMENTALANALYSIS_H
#define INCREMENTALANALYSIS_H

#include "NonLinearAnalysis.h"

class IncrementalAnalysis : public NonLinearAnalysis
{
protected: 
    int fNSteps=1;

public:
    IncrementalAnalysis() : NonLinearAnalysis(){};
    IncrementalAnalysis(CompMesh *cmesh, SolverType stype, double tol = 1.e-6, int maxIter = 10) : NonLinearAnalysis(cmesh,stype,tol,maxIter){};
    IncrementalAnalysis(Arlequin *arl, SolverType stype, double tol = 1.e-6, int maxIter = 10) : NonLinearAnalysis(arl,stype,tol,maxIter) {};
    int &NSteps(){return fNSteps;}

    void Run() override;
};

#endif