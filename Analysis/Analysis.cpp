#include "Analysis.h"
#include "Assemble.h"
#include "PETScSolver.h"
#include "EigenLinearSolver.h"
#include "GlobalMatrixT.h"
#include "PETScMatrix.h"

Analysis::Analysis(CompMesh *cmesh, SolverType stype){
    fMeshVector.resize(1);
    fMeshVector[0] = cmesh;
    fSolverType = stype;
    AllocateMonomodel();
};

Analysis::Analysis(Arlequin *arl, SolverType stype){
    fArlequin = arl;
    fMeshVector = fArlequin->MeshVec();
    fSolverType = stype;
    AllocateArlequin();        
};

Analysis::~Analysis() {
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

#ifdef HAS_PETSC
    fSolver = new PETScSolver(this);
#else
    fSolver = new EigenLinearSolver(this);
#endif
    fSolver->Solve();
}

void Analysis::AllocateMonomodel(){

    int numDOF = fMeshVector[0]->NGlobalDOF();
    std::cout << "Number of DOF = " << numDOF << std::endl;


#ifdef HAS_PETSC
    if (fSolverType == SolverType::EUmfpack){
        fGlobalMatrix = new GlobalMatrixT<PETScMatrix>(numDOF,numDOF);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, numDOF, numDOF, 100,NULL,&fGlobalStiffness);
    } else if (fSolverType == SolverType::ECholmod || fSolverType == SolverType::EKLU || fSolverType == SolverType::ESPQR){
        MatCreateSeqAIJ(PETSC_COMM_WORLD, numDOF, numDOF, 100,NULL,&fGlobalStiffness);
        MatSetOption(fGlobalStiffness, MAT_SYMMETRIC, PETSC_TRUE);
    } else {
        MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        numDOF, numDOF,100,NULL,300,NULL,&fGlobalStiffness); 
    }
    for (PetscInt i=0; i<numDOF; i++){
        double val = 1.e-20;
        MatSetValues(fGlobalStiffness,1,&i,1,&i,&val,ADD_VALUES);
    }
    //Create PETSc vectors
    VecCreate(PETSC_COMM_WORLD,&fGlobalRhs);
    VecSetSizes(fGlobalRhs,PETSC_DECIDE,numDOF);
    
    VecSetFromOptions(fGlobalRhs);
    VecDuplicate(fGlobalRhs,&fGlobalSolution);
#endif


}


void Analysis::AllocateArlequin(){

    int64_t numDOFGlobal = fMeshVector[0]->NGlobalDOF();
    int64_t numDOFLocal = fMeshVector[1]->NGlobalDOF();
    int64_t numDOFLagMul = fMeshVector[2]->NGlobalDOF();
    int64_t numDOF = numDOFGlobal + numDOFLocal + numDOFLagMul;
    if (fSolverType == SolverType::EUmfpack){
        MatCreateSeqAIJ(PETSC_COMM_WORLD, numDOF, numDOF, 150,NULL,&fGlobalStiffness);
    } else if (fSolverType == SolverType::ECholmod || fSolverType == SolverType::EKLU || fSolverType == SolverType::ESPQR){
        MatCreateSeqAIJ(PETSC_COMM_WORLD, numDOF, numDOF, 150,NULL,&fGlobalStiffness);
        MatSetOption(fGlobalStiffness, MAT_SYMMETRIC, PETSC_TRUE);
    } else {
        MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        numDOF, numDOF,100,NULL,300,NULL,&fGlobalStiffness); 
    }
    for (PetscInt i=0; i<numDOF; i++){
        double val = 1.e-20;
        MatSetValues(fGlobalStiffness,1,&i,1,&i,&val,ADD_VALUES);
    }

    //Create PETSc vectors
    VecCreate(PETSC_COMM_WORLD,&fGlobalRhs);
    VecSetSizes(fGlobalRhs,PETSC_DECIDE,numDOF);
    
    VecSetFromOptions(fGlobalRhs);
    VecDuplicate(fGlobalRhs,&fGlobalSolution);

    // MatAssemblyBegin(fGlobalStiffness,MAT_FINAL_ASSEMBLY);
    // MatAssemblyEnd(fGlobalStiffness,MAT_FINAL_ASSEMBLY);

    // VecAssemblyBegin(fGlobalRhs);
    // VecAssemblyEnd(fGlobalRhs);
}

void Analysis::PostProcessError(VecDouble &errorsTotal){

#ifdef HAS_PETSC
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
#endif

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

#ifdef HAS_PETSC
        if (errorsTotal.size()>0) MPI_Allreduce(&errorsProcess[0],&errorsTotal[0],errorsTotal.size(),MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
        if (rank == 0 && errorsTotal.size()>0){
            std::cout << "\n\nERROR REPORT - MESHVECTOR[" << imesh << "]:\n" << std::scientific << std::setprecision(10);
            for (int k = 0; k < errorsTotal.size(); k++){
                std::cout << "Errors[" << k << "] = " << sqrt(errorsTotal[k]) << "\n";
            }
        }
#else
        if (errorsTotal.size()>0) errorsTotal = errorsProcess;
        if (errorsTotal.size()>0){
            std::cout << "\n\nERROR REPORT - MESHVECTOR[" << imesh << "]:\n" << std::scientific << std::setprecision(10);
            for (int k = 0; k < errorsTotal.size(); k++){
                std::cout << "Errors[" << k << "] = " << sqrt(errorsTotal[k]) << "\n";
            }
        }
#endif
        
    }
}


