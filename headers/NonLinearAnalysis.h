#ifndef NONLINEARANALYSIS_H
#define NONLINEARANALYSIS_H

#include "LinearAnalysis.h"

class NonLinearAnalysis : public LinearAnalysis
{
private:
    /* data */
public:
    NonLinearAnalysis() : LinearAnalysis(){};
    NonLinearAnalysis(CompMesh *cmesh, SolverType stype) : LinearAnalysis(cmesh,stype){};
    NonLinearAnalysis(std::vector<CompMesh *> &meshvector, SolverType stype) : LinearAnalysis(meshvector,stype){};

    // void Compute() override;

    void UpdateSolution() override;

    void Run() override;
};

#endif