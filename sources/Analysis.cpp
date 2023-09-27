#include "Analysis.h"
#include "Element.h"
#include "Boundary.h"
#include "CompMesh.h"
#include "Assemble.h"

Analysis::~Analysis()
{
    KSPDestroy(&ksp); 
    VecDestroy(&fGlobalRhs); 
    VecDestroy(&fGlobalSolution); 
    MatDestroy(&fGlobalStiffness); 
}

void Analysis::Solve(){
    PetscErrorCode    ierr;
    //Create KSP context to solve the linear system
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
    
    ierr = KSPSetOperators(ksp,fGlobalStiffness,fGlobalStiffness);
    
    switch (fSolverType)
    {
    case SolverType::ESuiteSparse:
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCLU);
        PCFactorSetMatSolverType(pc, MATSOLVERUMFPACK);
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

    ierr = VecView(fGlobalSolution,PETSC_VIEWER_STDOUT_WORLD);
}

void Analysis::AllocateProblem(){

    PetscErrorCode    ierr;
    int numDOF = fMeshVector[0]->NumDOF();
    if (fSolverType == SolverType::ESuiteSparse){
        ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, numDOF, numDOF, 100,NULL,&fGlobalStiffness);
    } else {
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        numDOF, numDOF,100,NULL,300,NULL,&fGlobalStiffness); 
    }
    for (int i=0; i<numDOF; i++){
        double val = 1.e-20;
        ierr = MatSetValues(fGlobalStiffness,1,&i,1,&i,&val,ADD_VALUES);
    }

    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD,&fGlobalRhs);
    ierr = VecSetSizes(fGlobalRhs,PETSC_DECIDE,numDOF);
    
    ierr = VecSetFromOptions(fGlobalRhs);
    ierr = VecDuplicate(fGlobalRhs,&fGlobalSolution);
}


