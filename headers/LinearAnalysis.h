#ifndef LINEARANALYSIS_H
#define LINEARANALYSIS_H

#include "Analysis.h"

class LinearAnalysis : public Analysis {
public:
    LinearAnalysis(CompMesh *cmesh, SolverType stype) : Analysis(cmesh,stype){};
    ~LinearAnalysis(){
        delete [] this;
    };

    void Compute() override;

    void UpdateSolution() override;
};





#endif