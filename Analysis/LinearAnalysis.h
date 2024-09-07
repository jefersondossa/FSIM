#ifndef LINEARANALYSIS_HPP
#define LINEARANALYSIS_HPP

#include "Analysis.h"
#include "Assemble.h"

class LinearAnalysis : public Analysis {
public:
    LinearAnalysis() : Analysis(){};
    LinearAnalysis(CompMesh *cmesh, SolverType stype) : Analysis(cmesh,stype){};
    LinearAnalysis(Arlequin* arl, SolverType stype, bool reduced = true) : Analysis(arl,stype,reduced){};

    void Compute() override;
    void ComputeJacobian() override;
    void ComputeRhs() override;

    void UpdateSolution() override;

    ~LinearAnalysis(){};
};


#endif