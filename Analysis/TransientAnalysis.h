#ifndef TRANSIENTANALYSIS_H
#define TRANSIENTANALYSIS_H

#include "LinearAnalysis.h"
#include "NonLinearAnalysis.h"

class TransientAnalysis : public Analysis
{
public:
    TransientAnalysis() : Analysis(){};
    TransientAnalysis(CompMesh *cmesh, SolverType stype) : Analysis(cmesh,stype){};

    void Compute() override;

    void UpdateSolution() override;

    void Run() override;

};



#endif