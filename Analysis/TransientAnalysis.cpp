#include "TransientAnalysis.h"
#include "TransientWeakForm.h"
#include "VTUGenerator.h"

TransientAnalysis::TransientAnalysis(CompMesh *cmesh, SolverType stype, bool linear) : NonLinearAnalysis(cmesh,stype){
    IsLinear = linear;
};


void TransientAnalysis::Run(int64_t nsteps){

    for (int64_t i = 0; i < nsteps; i++){
        std::cout << "Time Step = " << i << std::endl;
        if (IsLinear){
            LinearAnalysis::Compute();
            LinearAnalysis::Solve();
            LinearAnalysis::UpdateSolution();
        } else {
            NonLinearAnalysis::Run();
        }
        // fGlobalMatrix->PrintMatrix();
        // fGlobalMatrix->PrintRhs();
        // fGlobalMatrix->PrintSolution();
        for (int i = 0; i < fMeshVector.size(); i++){
            // Loop over all weak forms in the mesh. If it is a transient weak form
            // update the time derivatives
            auto matvec = fMeshVector[i]->MaterialVector();
            for (auto const& [key, val] : matvec){
                TransientWeakForm *twf = dynamic_cast<TransientWeakForm*> (val);
                if(twf) twf->UpdateTimeDerivatives(fMeshVector[i]);
            }
        }
        fGlobalMatrix->ZeroMatrix();
        fGlobalMatrix->ZeroRhs();
        fGlobalMatrix->ZeroSolution();
        for (int64_t i=0; i<NEquations(); i++){
            double val = 1.e-20;
            fGlobalMatrix->AddValueMatrix(i,i,val);
        }
        VTUGenerator::PrintResults(this->MeshVector()[0],fFilename,fScalVars,fVectVars,i);
    }
}

void TransientAnalysis::PrintVariables(std::string filename, std::vector<std::string> &scal, std::vector<std::string> &vec){
    fScalVars = scal;
    fVectVars = vec;
    fFilename = filename;
}