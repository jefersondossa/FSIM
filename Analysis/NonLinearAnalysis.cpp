#include "NonLinearAnalysis.h"
#include "Assemble.h"

void NonLinearAnalysis::UpdateSolution(){
    VecScatter        ctx;
    Vec               All;
    PetscErrorCode    ierr;
    PetscScalar one = -1.;
    // VecAXPY(Solution(), one, fPreviousSolution);
    // VecCopy(Solution(), fPreviousSolution);

    ierr = VecDuplicate(Rhs(),&All);
    //Gathers the solution vector to the master process
    ierr = VecScatterCreateToAll(Solution(), &ctx, &All);
    ierr = VecScatterBegin(ctx, Solution(), All, INSERT_VALUES, SCATTER_FORWARD);
    ierr = VecScatterEnd(ctx, Solution(), All, INSERT_VALUES, SCATTER_FORWARD);
    ierr = VecScatterDestroy(&ctx);

    //Updates nodal values
    PetscInt Ione = 1;
    PetscInt Ii;
    PetscScalar val;
    
    int64_t nstartDOF = 0;
    for (int imesh = 0; imesh < this->MeshVector().size(); imesh++){
        if (imesh > 0) nstartDOF += this->MeshVector()[imesh-1]->NGlobalDOF();
        for (int i = 0; i < this->MeshVector()[imesh]->NNodes(); ++i){
            int nstate = this->MeshVector()[imesh]->NState();
            for (int k = 0; k<nstate; k++){
                Ii = nstartDOF + nstate*i+k;
                ierr = VecGetValues(All, Ione, &Ii, &val);
                this->MeshVector()[imesh]->NodeVec()[i] -> IncrementSolution(k,val);
            }
            // std::cout << "Solution mesh " << imesh << " " << i <<" " << val <<"\n";
        };
    }
    ierr = VecDestroy(&All); 
}

void NonLinearAnalysis::Run(){
    
    double NRL2norm = 1000.;
    int iteration = 0;

    VecDuplicate(this->Solution(),&fPreviousSolution);
    VecZeroEntries(fPreviousSolution);
    while (NRL2norm > fTolerance && iteration < fMaxIterations)
    {   
        std::clock_t t3 = std::clock();
        Compute();
        std::clock_t t4 = std::clock();
        std::cout << "Time assembling = " << 1000.*(t4-t3)/CLOCKS_PER_SEC/1000. << "s \n";
        // VecView(this->Rhs(),PETSC_VIEWER_STDOUT_WORLD);
        Solve();
        std::clock_t t5 = std::clock();
        std::cout << "Time Solving = " << 1000.*(t5-t4)/CLOCKS_PER_SEC/1000. << "s \n";
        // VecView(this->Solution(),PETSC_VIEWER_STDOUT_WORLD);
        // MatView(this->Stiffness(),PETSC_VIEWER_STDOUT_WORLD); 
        UpdateSolution();
        VecNorm(this->Solution(),NORM_2,&NRL2norm);
        std::cout << "Iteration " << iteration++ << ", Newton-Raphson residual = " << NRL2norm << std::endl;
        MatZeroEntries(this->Stiffness());
        VecZeroEntries(this->Rhs());
        VecZeroEntries(this->Solution());
        for (PetscInt i=0; i<NEquations(); i++){
            double val = 1.e-20;
            MatSetValues(this->Stiffness(),1,&i,1,&i,&val,ADD_VALUES);
        }
    }
    
    


}
