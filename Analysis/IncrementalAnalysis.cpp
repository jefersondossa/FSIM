#include "IncrementalAnalysis.h"
#include "VTUGenerator.h"

IncrementalAnalysis::IncrementalAnalysis(CompMesh *cmesh, SolverType stype, int nsteps, std::vector<L2Projection *> &bcinc, double tol, int maxIter) 
    : NonLinearAnalysis(cmesh,stype,tol,maxIter), fNSteps(nsteps){
    fIncrementBC = bcinc;
};


IncrementalAnalysis::IncrementalAnalysis(Arlequin *arl, SolverType stype, int nsteps, std::vector<L2Projection *> &bcinc, double tol, int maxIter) 
    : NonLinearAnalysis(arl,stype,tol,maxIter), fNSteps(nsteps){
    fIncrementBC = bcinc;
};

void IncrementalAnalysis::Run(){

    std::vector<VecDouble > increments(fIncrementBC.size());
    //Define the step size for each BC
    for (int i = 0; i < fIncrementBC.size(); i++){
        auto bcval = fIncrementBC[i]->BCValue();
        increments[i] = bcval/fNSteps;
    }
    
    int iStep = 0;
    while (iStep <= fNSteps){
        std::cout << "Computing step... " << iStep << std::endl;
        
        for (int i = 0; i < fIncrementBC.size(); i++){
            fIncrementBC[i]->BCValue() += increments[i] * iStep;
        }

        NonLinearAnalysis::Run();
        iStep++;
    }
}


void IncrementalAnalysis::Run(std::string filename, std::vector<std::string> &scalnames, std::vector<std::string> &vecnames){

    std::vector<VecDouble > increments(fIncrementBC.size());
    std::vector<VecDouble> initialbc(fIncrementBC.size());
    //Define the step size for each BC
    if (fIncrementValue.size() == 0){
        for (int i = 0; i < fIncrementBC.size(); i++){
            auto bcval = fIncrementBC[i]->BCValue();
            initialbc[i] = bcval;
            initialbc[i].setZero();
            increments[i] = bcval/fNSteps;
        }
    } else {
        for (int i = 0; i < fIncrementBC.size(); i++){
            initialbc[i] = fIncrementBC[i]->BCValue();
            increments[i] = fIncrementValue; 
        }
    }
    
    int iStep = 0;
    
    while (iStep <= fNSteps-1){
        std::cout << "Computing step... " << iStep << std::endl;
        
        if (fIncrementValue.size() == 0){
            for (int i = 0; i < fIncrementBC.size(); i++){
                fIncrementBC[i]->BCValue() = initialbc[i] + increments[i] * iStep;
            }
        } else {
            for (int i = 0; i < fIncrementBC.size(); i++){
                fIncrementBC[i]->BCValue() = initialbc[i] + increments[i] * iStep;
            }
        }

        NonLinearAnalysis::Run();
        VTUGenerator::PrintResults(MeshVector()[0],filename,scalnames,vecnames,iStep);
        iStep++;
    }
}
