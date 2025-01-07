#include "Analysis.h"
#include "Assemble.h"
#include "PETScSolver.h"
#include "EigenLinearSolver.h"
#include "ArlequinRedSolverEigen.h"
#include "ArlequinRedSolverPETSc.h"
#include "PETScMatrix.h"
#include "EigenSpMatrix.h"

Analysis::Analysis(CompMesh *cmesh, SolverType stype){
    fMeshVector.resize(1);
    fMeshVector[0] = cmesh;
    fReducedArlequin = false;
    fArlequin = nullptr;
    fSolverType = stype;
    AllocateMonomodel();
};

Analysis::Analysis(Arlequin *arl, SolverType stype, bool reduced){
    fReducedArlequin = reduced;
    fArlequin = arl;
    fMeshVector = fArlequin->MeshVec();
    fSolverType = stype;
    AllocateArlequin();
};


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
    if (fArlequin && fReducedArlequin){
        fSolver = new ArlequinRedSolverPETSc(this);
    } else {
        fSolver = new PETScSolver(this);
    }
#else
    if (fArlequin && fReducedArlequin){
        fSolver = new ArlequinRedSolverEigen(this);
    } else {
        fSolver = new EigenLinearSolver(this);
    }
#endif
    fSolver->Solve();
    delete fSolver;
    fSolver = nullptr;
}

void Analysis::AllocateMonomodel(){

    int numDOF = fMeshVector[0]->NGlobalDOF();
    std::cout << "Number of DOF = " << numDOF << std::endl;

    if(fGlobalMatrix)
    {
        delete fGlobalMatrix;
        fGlobalMatrix = nullptr;
    }

    if (fSolverType == SolverType::EUmfpack){
#ifdef HAS_PETSC
        fGlobalMatrix = new PETScMatrix(numDOF,numDOF,PETScMatType::ESeq);
#else 
        fGlobalMatrix = new EigenSpMatrix(numDOF,numDOF);
#endif
    } else if (fSolverType == SolverType::ECholmod || fSolverType == SolverType::EKLU || fSolverType == SolverType::ESPQR){
#ifdef HAS_PETSC
        fGlobalMatrix = new PETScMatrix(numDOF,numDOF,PETScMatType::ESeqSym);
#else
        fGlobalMatrix = new EigenSpMatrix(numDOF,numDOF);
#endif
    } else {
#ifdef HAS_PETSC
        fGlobalMatrix = new PETScMatrix(numDOF,numDOF,PETScMatType::EAij);
#else
        fGlobalMatrix = new EigenSpMatrix(numDOF,numDOF);
#endif
    }
    for (int64_t i=0; i<numDOF; i++){
        double val = 1.e-20;
        fGlobalMatrix->AddValueMatrix(i,i,val);
    }

}


void Analysis::AllocateArlequin(){
    if(fGlobalMatrix)
    {
        delete fGlobalMatrix;
        fGlobalMatrix = nullptr;
    }

    int64_t numDOFGlobal = fMeshVector[0]->NGlobalDOF();
    int64_t numDOFLocal = fMeshVector[1]->NGlobalDOF();
    int64_t numDOFLagMul = fMeshVector[2]->NGlobalDOF();
    int64_t numDOF = numDOFGlobal + numDOFLocal + numDOFLagMul;
    if (fSolverType == SolverType::EUmfpack){
#ifdef HAS_PETSC
    if (fReducedArlequin){
        fGlobalMatrix = new ArlequinMatRedPETSc(numDOFGlobal,numDOFLocal,numDOFLagMul,PETScMatType::ESeq);
    } else {
        fGlobalMatrix = new PETScMatrix(numDOF,numDOF,PETScMatType::ESeq);
    }
#else 
    if (fReducedArlequin){
        fGlobalMatrix = new ArlequinMatRedEigen(numDOFGlobal,numDOFLocal,numDOFLagMul);
    } else {
        fGlobalMatrix = new EigenSpMatrix(numDOF,numDOF);
    }
#endif
    } else if (fSolverType == SolverType::ECholmod || fSolverType == SolverType::EKLU || fSolverType == SolverType::ESPQR){
#ifdef HAS_PETSC
    if (fReducedArlequin){
        fGlobalMatrix = new ArlequinMatRedPETSc(numDOFGlobal,numDOFLocal,numDOFLagMul,PETScMatType::ESeqSym);
    } else {
        fGlobalMatrix = new PETScMatrix(numDOF,numDOF,PETScMatType::ESeqSym);
    }
#else
    if (fReducedArlequin){
        fGlobalMatrix = new ArlequinMatRedEigen(numDOFGlobal,numDOFLocal,numDOFLagMul);
    } else {
        fGlobalMatrix = new EigenSpMatrix(numDOF,numDOF);
    }
#endif
    } else {
#ifdef HAS_PETSC
    if (fReducedArlequin){
        fGlobalMatrix = new ArlequinMatRedPETSc(numDOFGlobal,numDOFLocal,numDOFLagMul,PETScMatType::EAij);
    } else {
        fGlobalMatrix = new PETScMatrix(numDOF,numDOF,PETScMatType::EAij);
    }
#else
    if (fReducedArlequin){
        fGlobalMatrix = new ArlequinMatRedEigen(numDOFGlobal,numDOFLocal,numDOFLagMul);
    } else {
        fGlobalMatrix = new EigenSpMatrix(numDOF,numDOF);
    }
#endif
    }
    for (int64_t i=0; i<numDOF; i++){
        double val = 1.e-20;
        fGlobalMatrix->AddValueMatrix(i,i,val);
    }
}

void Analysis::PostProcessError(VecDouble &errorsTotal){

#ifdef HAS_PETSC
    // PetscMPIInt rank;
    // PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
#endif
    int isArlequin = 0;
    if (fArlequin) isArlequin = 1;
    for (int imesh = 0; imesh < fMeshVector.size()-isArlequin; imesh++){
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
        // if (errorsTotal.size()>0) MPI_Allreduce(&errorsProcess[0],&errorsTotal[0],errorsTotal.size(),MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
        // if (rank == 0 && errorsTotal.size()>0){
        //     std::cout << "\n\nERROR REPORT - MESHVECTOR[" << imesh << "]:\n" << std::scientific << std::setprecision(10);
        //     for (int k = 0; k < errorsTotal.size(); k++){
        //         std::cout << "Errors[" << k << "] = " << sqrt(errorsTotal[k]) << "\n";
        //     }
        // }
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

void Analysis::Run(){
    std::cout << "Allocating problem..." << std::endl;
    std::clock_t t3 = std::clock();
    Compute();
    std::clock_t t4 = std::clock();
    std::cout << "Time assembling = " << 1000.*(t4-t3)/CLOCKS_PER_SEC/1000. << "s \n";
    Solve();
    std::clock_t t5 = std::clock();
    std::cout << "Time Solving = " << 1000.*(t5-t4)/CLOCKS_PER_SEC/1000. << "s \n";
    std::cout << "Updating solution..." << std::endl;
    UpdateSolution();
};
