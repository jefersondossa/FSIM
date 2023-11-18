#ifndef INCREMENTALANALYSIS_H
#define INCREMENTALANALYSIS_H

#include "NonLinearAnalysis.h"
#include "L2Projection.h"

class IncrementalAnalysis : public NonLinearAnalysis
{
protected: 
    int fNSteps = 1;
    std::vector<L2Projection *> fIncrementBC;

public:
    IncrementalAnalysis() : NonLinearAnalysis(){};
    IncrementalAnalysis(CompMesh *cmesh, SolverType stype, int nsteps, std::vector<L2Projection *> &bcinc, double tol = 1.e-6, int maxIter = 10);
    IncrementalAnalysis(Arlequin *arl, SolverType stype, int nsteps, std::vector<L2Projection *> &bcinc, double tol = 1.e-6, int maxIter = 10);

    int &NSteps(){return fNSteps;}

    void Run() override;
    void Run(std::string filename, std::vector<std::string> &scalnames, std::vector<std::string> &vecnames);
};

#endif