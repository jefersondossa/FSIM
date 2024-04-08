#include "LinearAnalysis.h"
#include "Assemble.h"

void LinearAnalysis::Compute(){
    if (this->MeshVector().size() == 1){
        Assemble::Monomodel(this);
    } else {
        Assemble::Arlequin(this);
    }
    this->GlobalMatrix()->MatVecAssemble();

    // this->GlobalMatrix()->PrintMatrix();
    // this->GlobalMatrix()->PrintRhs();
    // this->GlobalMatrix()->PrintSolution();
} 

void LinearAnalysis::ComputeJacobian(){
    if (this->MeshVector().size() == 1){
        Assemble::MonomodelMatrix(this);
    } else {
        Assemble::ArlequinMatrix(this);
    }
    this->GlobalMatrix()->MatAssemble();
    
    // this->GlobalMatrix()->PrintMatrix();
    // this->GlobalMatrix()->PrintRhs();
    // this->GlobalMatrix()->PrintSolution();
} 

void LinearAnalysis::ComputeRhs(){
    if (this->MeshVector().size() == 1){
        Assemble::MonomodelVector(this);
    } else {
        Assemble::ArlequinVector(this);
    }
    this->GlobalMatrix()->VecAssemble();
    
    // this->GlobalMatrix()->PrintMatrix();
    // this->GlobalMatrix()->PrintRhs();
    // this->GlobalMatrix()->PrintSolution();
} 

void LinearAnalysis::UpdateSolution(){
    this->GlobalMatrix()->ExpandSolution();

    //Updates nodal values
    int64_t Ii;
    double val;
    
    int64_t nstartDOF = 0;
    for (int imesh = 0; imesh < this->MeshVector().size(); imesh++){
        if (imesh > 0) nstartDOF += this->MeshVector()[imesh-1]->NGlobalDOF();
        for (int i = 0; i < this->MeshVector()[imesh]->NNodes(); ++i){
            int nstate = this->MeshVector()[imesh]->NodeVec()[i]->GetNStateVariables();
            for (int k = 0; k<nstate; k++){
                Ii = nstartDOF + nstate*i+k;
                val = this->GlobalMatrix()->GetValueSolution(Ii);
                // ierr = VecGetValues(All, Ione, &Ii, &val);
                this->MeshVector()[imesh]->NodeVec()[i] -> SetSolution(k,val);
            }
        };
    }
    this->GlobalMatrix()->ClearSolution();
}

