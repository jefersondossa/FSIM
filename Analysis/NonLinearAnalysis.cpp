#include "NonLinearAnalysis.h"
#include "Assemble.h"
#include "L2Projection.h"
#ifdef HAS_PETSC
#include "PETScSolver.h"
#endif

void NonLinearAnalysis::UpdateSolution(){
    
    this->GlobalMatrix()->ExpandSolution();

    //Updates nodal values
    int64_t Ione = 1;
    int64_t Ii;
    double val;
    
    int64_t nstartDOF = 0;
    for (int imesh = 0; imesh < this->MeshVector().size(); imesh++){
        if (imesh > 0) nstartDOF += this->MeshVector()[imesh-1]->NGlobalDOF();
        for (int i = 0; i < this->MeshVector()[imesh]->NNodes(); ++i){
            int nstate = this->MeshVector()[imesh]->NState();
            for (int k = 0; k<nstate; k++){
                Ii = nstartDOF + nstate*i+k;
                val = this->GlobalMatrix()->GetValueSolution(Ii);
                // ierr = VecGetValues(All, Ione, &Ii, &val);
                this->MeshVector()[imesh]->NodeVec()[i] -> IncrementSolution(k,val);
            }
            // std::cout << "Solution mesh " << imesh << " " << i <<" " << val <<"\n";
        };
    }
    this->GlobalMatrix()->ClearSolution();
}

#ifdef HAS_PETSC

PetscErrorCode NonLinearAnalysis::FormFunction(SNES snes, Vec u, Vec b, void *ptr){
    NonLinearAnalysis *an = static_cast<NonLinearAnalysis * >(ptr);
    an->GlobalMatrix()->ExpandSolution();

    //Updates nodal values
    int64_t Ione = 1;
    int64_t Ii;


    // std::cout << "PrintSol = \n"; 
    // VecView(u,PETSC_VIEWER_STDOUT_WORLD);
    // an->GlobalMatrix()->PrintSolution();
    

    VecScatter  ctx;
    Vec         SolAll;
    //Gathers the solution vector to the master process
    VecScatterCreateToAll(u, &ctx, &SolAll);
    VecScatterBegin(ctx, u, SolAll, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, u, SolAll, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&ctx);
    // std::cout << "PrintSolAll = \n"; 
    // VecView(SolAll,PETSC_VIEWER_STDOUT_WORLD);
    
    int64_t nstartDOF = 0;
    for (int imesh = 0; imesh < an->MeshVector().size(); imesh++){
        if (imesh > 0) nstartDOF += an->MeshVector()[imesh-1]->NGlobalDOF();
        for (int i = 0; i < an->MeshVector()[imesh]->NNodes(); ++i){
            int nstate = an->MeshVector()[imesh]->NState();
            for (int k = 0; k<nstate; k++){
                Ii = nstartDOF + nstate*i+k;
                int64_t Ione = 1;
                PetscScalar val;
                VecGetValues(SolAll, Ione, &Ii, &val);
                // double prevsol = an->MeshVector()[imesh]->NodeVec()[i] -> GetSolution(k);
                an->MeshVector()[imesh]->NodeVec()[i] -> SetSolution(k,val);
            }
        };
    }

    VecDestroy(&SolAll);
    // this->GlobalMatrix()->ClearSolution();
    an->GlobalMatrix()->ZeroRhs();
    an->ComputeRhs();
    // std::cout << "FORM FUNCTION Before \n "<< std::endl;
    // VecView(b,PETSC_VIEWER_STDOUT_WORLD);
    
    PETScMatrix *petscmat = dynamic_cast<PETScMatrix *> (an->GlobalMatrix());
    PetscScalar mone = -1.;
    VecScale(petscmat->Rhs(),mone);
    // VecCopy(petscmat->Rhs(),b);
    // VecCopy(petscmat->Rhs(),b);
    VecCopy(u,petscmat->Solution());
    // std::cout << "FORM FUNCTION RHS \n "<< std::endl;
    // an->GlobalMatrix()->PrintRhs();
    // std::cout << "FORM FUNCTION B \n "<< std::endl;
    // VecView(b,PETSC_VIEWER_STDOUT_WORLD);
    

    an->GlobalMatrix()->ClearSolution();
    return 0;
}


PetscErrorCode NonLinearAnalysis::FormJacobian(SNES snes,Vec u,Mat A, Mat B,void *ptr){
    NonLinearAnalysis *an = static_cast<NonLinearAnalysis * >(ptr);
    // an->GlobalMatrix()->ZeroMatrix();
    an->ComputeJacobian();
    return 0;
}

#endif

void NonLinearAnalysis::Run(){
    
    double NRL2norm = 1000.;
    int iteration = 0;
    
    // //Save Previous solution
    // for (int imesh = 0; imesh < this->MeshVector().size(); imesh++){
    //     for (int64_t inode = 0; inode < this->MeshVector()[imesh]->NNodes(); inode++){
    //         int nstate = this->MeshVector()[imesh]->NState();
    //         for (int istate = 0; istate < nstate; istate++){
    //             double Sol = this->MeshVector()[imesh]->NodeVec()[inode]->GetSolution(istate);
    //             this->MeshVector()[imesh]->NodeVec()[inode]->SetPreviousSolution(istate,Sol);
    //         }
    //     }
    // }
#ifdef HAS_PETSC
    
    PETScMatrix *petscmat = dynamic_cast<PETScMatrix *> (this->GlobalMatrix());
    fSolver = new PETScSolver(this);
    PETScSolver *petscsol = dynamic_cast<PETScSolver *> (this->Solver());

    SNESCreate(PETSC_COMM_WORLD,&fSNES);
    
    SNESSetFunction(fSNES,petscmat->Rhs(),FormFunction,&(*this));
    SNESSetJacobian(fSNES,petscmat->Matrix(),petscmat->Matrix(),FormJacobian,&(*this));

    SNESGetKSP(fSNES,&(*(petscsol->KSPSolver())));
    
    SNESSetFromOptions(fSNES);        
    // SNESView(fSNES,PETSC_VIEWER_STDOUT_WORLD);

    PetscInt iterations = 0;
    PetscReal fgnorm = 0.;
    PetscViewerAndFormat *vf;

    SNESSetForceIteration(fSNES,PETSC_TRUE);
    SNESSetType(fSNES, SNESNEWTONLS);
    SNESSolve(fSNES,NULL,petscmat->Solution());
    // petscmat->PrintSolution();
    // SNESMonitorSet(fSNES,Monitor,NULL,NULL);
    SNESGetIterationNumber(fSNES,&iterations);
    SNESDestroy(&fSNES);

    UpdateSolution();

    this->GlobalMatrix()->ZeroMatrix();
    this->GlobalMatrix()->ZeroRhs();
    // this->GlobalMatrix()->ZeroSolution();

#else
    //Iterative Process
    while (NRL2norm > fTolerance && iteration < fMaxIterations)
    {   
        std::ofstream output("plasticity.txt",std::ios::app);
        output << "ITERATION = " << iteration << std::endl;
        std::clock_t t3 = std::clock();
        Compute();
        std::clock_t t4 = std::clock();
        std::cout << "Time assembling = " << 1000.*(t4-t3)/CLOCKS_PER_SEC/1000. << "s \n";
        Solve();
        std::clock_t t5 = std::clock();
        std::cout << "Time Solving = " << 1000.*(t5-t4)/CLOCKS_PER_SEC/1000. << "s \n";
        UpdateSolution();
        // this->GlobalMatrix()->PrintMatrix();
        this->GlobalMatrix()->PrintRhs();
        this->GlobalMatrix()->PrintSolution();
        NRL2norm = this->GlobalMatrix()->SolutionNorm();
        std::cout << "Iteration " << iteration++ << ", Newton-Raphson residual = " << NRL2norm << std::endl;
        this->GlobalMatrix()->ZeroMatrix();
        this->GlobalMatrix()->ZeroRhs();
        this->GlobalMatrix()->ZeroSolution();
        for (int64_t i=0; i<NEquations(); i++){
            double val = 1.e-20;
            this->GlobalMatrix()->AddValueMatrix(i,i,val);
        }
    }
#endif

}
