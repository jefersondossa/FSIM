#ifndef NONLINEARANALYSIS_H
#define NONLINEARANALYSIS_H

#include "Analysis.h"
// class CompMesh;

class NonLinearAnalysis : public Analysis
{
private:
    /* data */
public:
    NonLinearAnalysis(CompMesh *cmesh, SolverType stype):Analysis(cmesh,stype){};
    ~NonLinearAnalysis();
};


NonLinearAnalysis::~NonLinearAnalysis()
{
}



#endif