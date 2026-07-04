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
    int Ii;
    REAL val;
    
    

    if (fMeshVector[0]->NConnects() == 0){
        for (auto cel : this->MeshVector()[0]->ElementVec()){
            if (!cel) continue;
            std::vector<Connect *> c = cel->getConnectivity();
            for (int  iconnect = 0; iconnect < c.size(); iconnect++){
                Connect * con = c[iconnect];
                int nstate = con->GetNStateVariables();
                for (int k = 0; k<nstate; k++){
                    Ii = con->GetSequenceNumber() + k;
                    if (Ii < 0) continue;
                    val = this->GlobalMatrix()->GetValueSolution(Ii);
                    // ierr = VecGetValues(All, Ione, &Ii, &val);
                    int nshape = con->GetNShapeFunctions();
                    if (nshape == 0) continue;
                    REAL sol = con -> Solution()[k];
                    con -> SetPreviousSolution(k,sol);
                    con -> SetSolution(k,val);
                }
            }
        }
    } else {
        int nstartDOF = 0;
        for (int imesh = 0; imesh < this->MeshVector().size(); imesh++){
            if (imesh > 0) nstartDOF += this->MeshVector()[imesh-1]->NGlobalDOF();
            
            for (int iconnect = 0; iconnect < this->MeshVector()[imesh]->NConnects(); iconnect++){
                Connect * c = this->MeshVector()[imesh]->ConnectVec()[iconnect];
                int fSeqnum = c->GetSequenceNumber();
                if (fSeqnum < 0) continue;
                int nstate = c->GetNStateVariables();
                if (c->GetNShapeFunctions() == 0) continue;
                for (int k = 0; k<nstate; k++){
                    Ii = nstartDOF + fSeqnum + k;
                    val = this->GlobalMatrix()->GetValueSolution(Ii);
                    // ierr = VecGetValues(All, Ione, &Ii, &val);
                    int nshape = c->GetNShapeFunctions();
                    if (nshape == 0) continue;
                    REAL sol = c -> Solution()[k];
                    c -> SetPreviousSolution(k,sol);
                    c -> SetSolution(k,val);
                }
            };
        }
    }

    // for (int iconnect = 0; iconnect < fMeshVector[0]->NConnects(); iconnect++){
    //     Connect * c = fMeshVector[0]->ConnectVec()[iconnect];
    //     int fSeqnum = c->GetSequenceNumber();
    //     if (fSeqnum < 0) continue;
    //     int nstate = c->GetNStateVariables();
    //     for (int k = 0; k<nstate; k++){
    //         Ii = fSeqnum + k;
    //         val = this->GlobalMatrix()->GetValueSolution(Ii);
    //         // ierr = VecGetValues(All, Ione, &Ii, &val);
    //         int nshape = c->GetNShapeFunctions();
    //         if (nshape == 0) continue;
    //         REAL sol = c -> Solution()[k];
    //         c -> SetPreviousSolution(k,sol);
    //         c -> SetSolution(k,val);
    //     }
    // }


    this->GlobalMatrix()->ClearSolution();
}

void LinearAnalysis::Run(){
    std::cout << "Allocating problem..." << std::endl;
    std::clock_t t3 = std::clock();
    Compute();
    std::clock_t t4 = std::clock();
    std::cout << "Time assembling = " << 1000.*(t4-t3)/CLOCKS_PER_SEC/1000. << "s \n";
    // fGlobalMatrix->PrintMatrix();
    // fGlobalMatrix->PrintRhs();
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

