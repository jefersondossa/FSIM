#ifndef TRANSIENTANALYSIS_H
#define TRANSIENTANALYSIS_H

#include "LinearAnalysis.h"
#include "NonLinearAnalysis.h"

class TransientAnalysis : public Analysis
{
private:
    bool IsLinear;

public:
    TransientAnalysis() : Analysis(){};
    TransientAnalysis(CompMesh *cmesh, SolverType stype, bool linear = true);

    void Compute() override;

    void UpdateSolution() override;

    void Run() override;

};



#endif