#include "TransientAnalysis.h"
#include "TransientWeakForm.h"
#include "VTUGenerator.h"
#include "L2Projection.h"

TransientAnalysis::TransientAnalysis(CompMesh *cmesh, SolverType stype, bool linear) : NonLinearAnalysis(cmesh,stype){
    IsLinear = linear;
};


void TransientAnalysis::Run(int64_t nsteps){

    for (int64_t i = 0; i < nsteps; i++){
        
        std::cout << "Time Step = " << i << std::endl;
        if (IsLinear){
            fGlobalMatrix->ZeroMatrix();
            fGlobalMatrix->ZeroRhs();
            fGlobalMatrix->ZeroSolution();
            LinearAnalysis::Run();
        } else {
            NonLinearAnalysis::Run();
        }
        for (int imesh = 0; imesh < fMeshVector.size(); imesh++){
            // Loop over all weak forms in the mesh. If it is a transient weak form
            // update the time derivatives
            auto matvec = fMeshVector[imesh]->MaterialVector();
            for (auto const& [key, val] : matvec){
                L2Projection *l2p = dynamic_cast<L2Projection*> (val);
                if (l2p){
                    l2p->SetTimeInstant(i+1);
                    continue; // Skip L2Projection, it does not have time derivatives
                } 
                TransientWeakForm *twf = dynamic_cast<TransientWeakForm*> (val);
                if(twf) {
                    twf->UpdateTimeDerivatives(fMeshVector[imesh]);
                    twf->SetTimeInstant(i+1);
                };
            }
        }
        
        VTUGenerator::PrintResults(this->MeshVector()[0],fFilename,fScalVars,fVectVars,{},i);
        fGlobalMatrix->ZeroMatrix();
        fGlobalMatrix->ZeroRhs();
        fGlobalMatrix->ZeroSolution();
        for (int64_t i=0; i<NEquations(); i++){
            double val = 1.e-20;
            fGlobalMatrix->AddValueMatrix(i,i,val);
        }
    }
}

void TransientAnalysis::PrintVariables(std::string filename, std::vector<std::string> &scal, std::vector<std::string> &vec){
    fScalVars = scal;
    fVectVars = vec;
    fFilename = filename;
}