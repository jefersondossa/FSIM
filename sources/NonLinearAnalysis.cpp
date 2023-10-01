#include "NonLinearAnalysis.h"
#include "Assemble.h"

// void NonLinearAnalysis::Compute(){
//     if (this->MeshVector().size() == 1){
//         Assemble::Monomodel(this);
//     } else {
//         Assemble::Arlequin(this);
//     }
//     PetscErrorCode    ierr;

//     //Assemble matrices and vectors
//     ierr = MatAssemblyBegin(this->Stiffness(),MAT_FINAL_ASSEMBLY);
//     ierr = MatAssemblyEnd(this->Stiffness(),MAT_FINAL_ASSEMBLY);
    
//     ierr = VecAssemblyBegin(this->Rhs());
//     ierr = VecAssemblyEnd(this->Rhs());

//     MatView(this->Stiffness(),PETSC_VIEWER_STDOUT_WORLD);
//     MatView(this->Stiffness(),PETSC_VIEWER_DRAW_WORLD);
//     VecView(this->Rhs(),PETSC_VIEWER_STDOUT_WORLD);
// } 

void NonLinearAnalysis::UpdateSolution(){
    VecScatter        ctx;
    Vec               All, Allu;
    PetscErrorCode    ierr;

    ierr = VecDuplicate(Rhs(),&All);
    //Gathers the solution vector to the master process
    ierr = VecScatterCreateToAll(Solution(), &ctx, &All);
    ierr = VecScatterBegin(ctx, Solution(), All, INSERT_VALUES, SCATTER_FORWARD);
    ierr = VecScatterEnd(ctx, Solution(), All, INSERT_VALUES, SCATTER_FORWARD);
    ierr = VecScatterDestroy(&ctx);

    ierr = VecScatterCreateToAll(Rhs(), &ctx, &Allu);
    ierr = VecScatterBegin(ctx, Rhs(), Allu, INSERT_VALUES, SCATTER_FORWARD);
    ierr = VecScatterEnd(ctx, Rhs(), Allu, INSERT_VALUES, SCATTER_FORWARD);
    ierr = VecScatterDestroy(&ctx);
    
    //Updates nodal values
    PetscInt Ione = 1;
    PetscInt Ii;
    PetscScalar val;
    
    int64_t nstartDOF = 0;
    for (int imesh = 0; imesh < this->MeshVector().size(); imesh++){
        if (imesh > 0) nstartDOF += this->MeshVector()[imesh-1]->NGlobalDOF();
        for (int i = 0; i < this->MeshVector()[imesh]->NNodes(); ++i){
            int nstate = this->MeshVector()[imesh]->NodeVec()[i]->GetNStateVariables();
            for (int k = 0; k<nstate; k++){
                Ii = nstartDOF + nstate*i+k;
                ierr = VecGetValues(All, Ione, &Ii, &val);
                this->MeshVector()[imesh]->NodeVec()[i] -> IncrementSolution(k,val);
            }
            // std::cout << "Solution mesh " << imesh << " " << i <<" " << val <<"\n";
        };
    }
    ierr = VecDestroy(&All); 
    ierr = VecDestroy(&Allu); 
}

void NonLinearAnalysis::Run(){
    
    double NRL2norm = 1000.;
    double NRtolerance = 1.e-6;
    int iteration = 0;

    while (NRL2norm > NRtolerance)
    {   
        Compute();
        Solve();
        // VecView(this->Solution(),PETSC_VIEWER_STDOUT_WORLD);
        UpdateSolution();
        VecNorm(this->Solution(),NORM_2,&NRL2norm);
        std::cout << "Iteration " << iteration++ << ", Newton-Raphson residual = " << NRL2norm << std::endl;
        // if (iteration == 2) break;
    }
    
    


}
