#include "TransientAnalysis.h"

TransientAnalysis::TransientAnalysis(CompMesh *cmesh, SolverType stype, bool linear) : NonLinearAnalysis(cmesh,stype){
    IsLinear = linear;
};


void TransientAnalysis::Run(int64_t nsteps){

    for (int64_t i = 0; i < nsteps; i++){
        if (IsLinear){
            LinearAnalysis::Run();
        } else {
            NonLinearAnalysis::Run();
        } 
    }
}