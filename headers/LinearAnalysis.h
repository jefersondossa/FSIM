#ifndef LINEARANALYSIS_HPP
#define LINEARANALYSIS_HPP

#include "Analysis.h"
#include "Assemble.h"

class LinearAnalysis : public Analysis {
public:
    LinearAnalysis() : Analysis(){};
    LinearAnalysis(CompMesh *cmesh, SolverType stype) : Analysis(cmesh,stype){};
    LinearAnalysis(std::vector<CompMesh *> &meshvector, SolverType stype) : Analysis(meshvector,stype){};

    void Compute() override;

    void UpdateSolution() override;

    ~LinearAnalysis(){};
};


#endif