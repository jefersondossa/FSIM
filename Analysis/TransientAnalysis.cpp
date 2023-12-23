#include "TransientAnalysis.h"

TransientAnalysis::TransientAnalysis(CompMesh *cmesh, SolverType stype, bool linear) : Analysis(cmesh,stype){
    IsLinear = linear;
};
