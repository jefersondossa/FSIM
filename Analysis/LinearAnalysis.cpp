#include "LinearAnalysis.h"
#include "Assemble.h"

void LinearAnalysis::Compute(){
    if (this->MeshVector().size() == 1){
        Assemble::Monomodel(this);
    } else {
        Assemble::Arlequin(this);
    }
    this->GlobalMatrix()->MatVecAssemble();
} 

void LinearAnalysis::ComputeJacobian(){
    if (this->MeshVector().size() == 1){
        Assemble::MonomodelMatrix(this);
    } else {
        Assemble::ArlequinMatrix(this);
    }
    this->GlobalMatrix()->MatAssemble();
} 

void LinearAnalysis::ComputeRhs(){
    if (this->MeshVector().size() == 1){
        Assemble::MonomodelVector(this);
    } else {
        Assemble::ArlequinVector(this);
    }
    this->GlobalMatrix()->VecAssemble();
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
            int nstate = this->MeshVector()[imesh]->ConnectVec()[i]->GetNStateVariables();
            for (int k = 0; k<nstate; k++){
                Ii = nstartDOF + nstate*i+k;
                val = this->GlobalMatrix()->GetValueSolution(Ii);
                // ierr = VecGetValues(All, Ione, &Ii, &val);
                double sol = this->MeshVector()[imesh]->ConnectVec()[i] -> Solution()[k];
                this->MeshVector()[imesh]->ConnectVec()[i] -> SetPreviousSolution(k,sol);
                this->MeshVector()[imesh]->ConnectVec()[i] -> SetSolution(k,val);
            }
        };
    }
    this->GlobalMatrix()->ClearSolution();
}

void LinearAnalysis::Run(){
    std::cout << "Allocating problem..." << std::endl;
    std::clock_t t3 = std::clock();
    Compute();
    std::clock_t t4 = std::clock();
    std::cout << "Time assembling = " << 1000.*(t4-t3)/CLOCKS_PER_SEC/1000. << "s \n";
    fGlobalMatrix->PrintMatrix();
    fGlobalMatrix->PrintRhs();
    // fGlobalMatrix->PrintSolution();
    Solve();
    std::clock_t t5 = std::clock();
    std::cout << "Time Solving = " << 1000.*(t5-t4)/CLOCKS_PER_SEC/1000. << "s \n";
    std::cout << "Updating solution..." << std::endl;
    LinearAnalysis::UpdateSolution();
    // fGlobalMatrix->PrintMatrix();
    // fGlobalMatrix->PrintRhs();
    // fGlobalMatrix->PrintSolution();
};

