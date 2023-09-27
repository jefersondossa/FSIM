#ifndef NONLINEARANALYSIS_H
#define NONLINEARANALYSIS_H

#include "Analysis.h"
// class CompMesh;

class NonLinearAnalysis : public Analysis
{
private:
    /* data */
public:
    NonLinearAnalysis() : Analysis(){};
    NonLinearAnalysis(CompMesh *cmesh, SolverType stype) : Analysis(cmesh,stype){};

    void Compute() override;

    void UpdateSolution() override;

    void Run() override;

    ~NonLinearAnalysis();
};


NonLinearAnalysis::~NonLinearAnalysis()
{
}



#endif