#ifndef TRANSIENTANALYSIS_H
#define TRANSIENTANALYSIS_H

#include "LinearAnalysis.h"
#include "NonLinearAnalysis.h"

class TransientAnalysis : public NonLinearAnalysis
{
private:
    bool IsLinear;

public:
    TransientAnalysis() : NonLinearAnalysis(){};
    TransientAnalysis(CompMesh *cmesh, SolverType stype, bool linear = true);

    void Run(int64_t nsteps);

    void Run() override{
        PanicButton();
    }

};



#endif