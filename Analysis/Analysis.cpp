#include "Analysis.h"
#include "Assemble.h"

Analysis::Analysis(Arlequin *arl, SolverType stype){
    fArlequin = arl;
    fMeshVector = fArlequin->MeshVec();
    fSolverType = stype;
    AllocateArlequin();        
};

Analysis::~Analysis()
{
    KSPDestroy(&ksp); 
    VecDestroy(&fGlobalRhs); 
    VecDestroy(&fGlobalSolution); 
    MatDestroy(&fGlobalStiffness); 
}

int64_t Analysis::NEquations(){
    int64_t fNumEquations = 0;
    for (auto mesh:fMeshVector){
        fNumEquations += mesh->NGlobalDOF();
    }
    
    return fNumEquations;
}

void Analysis::Solve(){
    std::cout << "Solving..." << std::endl;
    PetscErrorCode    ierr;
    //Create KSP context to solve the linear system
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
    
    ierr = KSPSetOperators(ksp,fGlobalStiffness,fGlobalStiffness);
    
    switch (fSolverType)
    {
    case SolverType::EUmfpack:
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCLU);
        PCFactorSetMatSolverType(pc, MATSOLVERUMFPACK);
        break;
    case SolverType::ECholmod:
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCCHOLESKY);
        PCFactorSetMatSolverType(pc, MATSOLVERCHOLMOD);
        break;
    case SolverType::EMumps:
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCLU);
        PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
        break;

    case SolverType::EIterative:
        KSPSetType(ksp,KSPFGMRES);
        KSPGetPC(ksp, &pc);
        PCSetType(pc,PCBJACOBI);
        KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,200);
        break;

    default:
        PanicButton();
        break;
    }

    //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    
    ierr = KSPSolve(ksp,fGlobalRhs,fGlobalSolution);
    // ierr = KSPGetTotalIterations(ksp, &iterations); 

    // ierr = VecView(fGlobalSolution,PETSC_VIEWER_STDOUT_WORLD);
}

void Analysis::AllocateMonomodel(){

    PetscErrorCode    ierr;
    int numDOF = fMeshVector[0]->NGlobalDOF();
    if (fSolverType == SolverType::EUmfpack || fSolverType == SolverType::ECholmod){
        ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, numDOF, numDOF, 100,NULL,&fGlobalStiffness);
    } else {
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        numDOF, numDOF,100,NULL,300,NULL,&fGlobalStiffness); 
    }
    for (PetscInt i=0; i<numDOF; i++){
        double val = 1.e-20;
        ierr = MatSetValues(fGlobalStiffness,1,&i,1,&i,&val,ADD_VALUES);
    }

    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD,&fGlobalRhs);
    ierr = VecSetSizes(fGlobalRhs,PETSC_DECIDE,numDOF);
    
    ierr = VecSetFromOptions(fGlobalRhs);
    ierr = VecDuplicate(fGlobalRhs,&fGlobalSolution);
}


void Analysis::AllocateArlequin(){

    PetscErrorCode    ierr;
    int64_t numDOFGlobal = fMeshVector[0]->NGlobalDOF();
    int64_t numDOFLocal = fMeshVector[1]->NGlobalDOF();
    int64_t numDOFLagMul = fMeshVector[2]->NGlobalDOF();
    int64_t numDOF = numDOFGlobal + numDOFLocal + numDOFLagMul;
    if (fSolverType == SolverType::EUmfpack || fSolverType == SolverType::ECholmod){
        ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, numDOF, numDOF, 150,NULL,&fGlobalStiffness);
    } else {
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        numDOF, numDOF,100,NULL,300,NULL,&fGlobalStiffness); 
    }
    for (PetscInt i=0; i<numDOF; i++){
        double val = 1.e-20;
        ierr = MatSetValues(fGlobalStiffness,1,&i,1,&i,&val,ADD_VALUES);
    }

    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD,&fGlobalRhs);
    ierr = VecSetSizes(fGlobalRhs,PETSC_DECIDE,numDOF);
    
    ierr = VecSetFromOptions(fGlobalRhs);
    ierr = VecDuplicate(fGlobalRhs,&fGlobalSolution);

    // ierr = MatAssemblyBegin(fGlobalStiffness,MAT_FINAL_ASSEMBLY);
    // ierr = MatAssemblyEnd(fGlobalStiffness,MAT_FINAL_ASSEMBLY);

    // ierr = VecAssemblyBegin(fGlobalRhs);
    // ierr = VecAssemblyEnd(fGlobalRhs);
}

void Analysis::PostProcessError(VecDouble &errorsTotal){

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    for (int imesh = 0; imesh < fMeshVector.size(); imesh++){
        VecDouble errorsProcess;
        // Loop over the elements
        for (int jel = fMeshVector[imesh]->NElements(); jel--; ){

            VecDouble errors;

            fMeshVector[imesh]->ElementVec()[jel] -> ComputeError(errors);
            if (errors.size() == 0) continue;
            errorsProcess.resize(errors.size());
            errorsProcess += errors;

        }; //Elements
        errorsTotal.resize(errorsProcess.size());
        errorsTotal.setZero();

        if (errorsTotal.size()>0) MPI_Allreduce(&errorsProcess[0],&errorsTotal[0],errorsTotal.size(),MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

        if (rank == 0 && errorsTotal.size()>0){
            std::cout << "\n\nERROR REPORT - MESHVECTOR[" << imesh << "]:\n" << std::scientific << std::setprecision(10);
            for (int k = 0; k < errorsTotal.size(); k++){
                std::cout << "Errors[" << k << "] = " << sqrt(errorsTotal[k]) << "\n";
            }
        }
    }
}


