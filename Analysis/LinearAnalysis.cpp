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
    
    // int64_t nstartDOF = 0;
    // for (int imesh = 0; imesh < this->MeshVector().size(); imesh++){
    //     if (imesh > 0) nstartDOF += this->MeshVector()[imesh-1]->NGlobalDOF();
    //     for (int i = 0; i < this->MeshVector()[imesh]->NConnects(); ++i){
    //         for (int iconnect = 0; iconnect < this->MeshVector()[imesh]->NConnects(); iconnect++){
    //             int64_t fSeqnum = this->MeshVector()[imesh]->ConnectVec()[iconnect]->GetSequenceNumber();
    //             if (fSeqnum < 0) continue;
    //             int nstate = this->MeshVector()[imesh]->ConnectVec()[iconnect]->GetNStateVariables();
    //             for (int k = 0; k<nstate; k++){
    //                 Ii = nstartDOF + fSeqnum + k;
    //                 val = this->GlobalMatrix()->GetValueSolution(Ii);
    //                 // ierr = VecGetValues(All, Ione, &Ii, &val);
    //                 int nshape = this->MeshVector()[imesh]->ConnectVec()[iconnect]->GetNShapeFunctions();
    //                 if (nshape == 0) continue;
    //                 double sol = this->MeshVector()[imesh]->ConnectVec()[iconnect] -> Solution()[k];
    //                 this->MeshVector()[imesh]->ConnectVec()[iconnect] -> SetPreviousSolution(k,sol);
    //                 this->MeshVector()[imesh]->ConnectVec()[iconnect] -> SetSolution(k,val);
    //             }
                
    //         }
    //     };
    // }

    for (int iconnect = 0; iconnect < fMeshVector[0]->NConnects(); iconnect++){
        Connect * c = fMeshVector[0]->ConnectVec()[iconnect];
        int64_t fSeqnum = c->GetSequenceNumber();
        if (fSeqnum < 0) continue;
        int nstate = c->GetNStateVariables();
        for (int k = 0; k<nstate; k++){
            Ii = fSeqnum + k;
            val = this->GlobalMatrix()->GetValueSolution(Ii);
            // ierr = VecGetValues(All, Ione, &Ii, &val);
            int nshape = c->GetNShapeFunctions();
            if (nshape == 0) continue;
            double sol = c -> Solution()[k];
            c -> SetPreviousSolution(k,sol);
            c -> SetSolution(k,val);
        }
        
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
    fGlobalMatrix->PrintSolution();
    Solve();
    std::clock_t t5 = std::clock();
    std::cout << "Time Solving = " << 1000.*(t5-t4)/CLOCKS_PER_SEC/1000. << "s \n";
    std::cout << "Updating solution..." << std::endl;
    LinearAnalysis::UpdateSolution();
    // fGlobalMatrix->PrintMatrix();
    // fGlobalMatrix->PrintRhs();
    // fGlobalMatrix->PrintSolution();
};

