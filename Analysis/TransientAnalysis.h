#ifndef TRANSIENTANALYSIS_H
#define TRANSIENTANALYSIS_H

#include "LinearAnalysis.h"
#include "NonLinearAnalysis.h"

class TransientAnalysis : public NonLinearAnalysis
{
private:
    bool IsLinear;
    std::vector<std::string> fScalVars;
    std::vector<std::string> fVectVars;
    std::string fFilename;
    int fNSteps;

public:
    TransientAnalysis() : NonLinearAnalysis(){};
    TransientAnalysis(CompMesh *cmesh, SolverType stype, bool linear = true);

    void Run(int nsteps);

    void Run() override{
        PanicButton();
    }
    void PrintVariables(std::string filename, std::vector<std::string> &scal, std::vector<std::string> &vec);

};



#endif