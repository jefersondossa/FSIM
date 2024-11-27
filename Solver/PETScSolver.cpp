#include "PETScSolver.h"

#ifdef HAS_PETSC

PETScSolver::PETScSolver(Analysis *an) : LinearSolver(an){
    //Create KSP context to solve the linear system
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    auto * pmat = dynamic_cast<PETScMatrix*> (fAnalysis->GlobalMatrix());
    KSPSetOperators(ksp,pmat->Matrix(),pmat->Matrix());
}

PETScSolver::~PETScSolver(){
    KSPDestroy(&ksp);
}

void PETScSolver::Solve(){

    switch (fAnalysis->SType())
    {
    case SolverType::EUmfpack:
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCLU);
        PCFactorSetMatSolverType(pc, MATSOLVERUMFPACK);
        break;
    case SolverType::EKLU:
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCLU);
        PCFactorSetMatSolverType(pc, MATSOLVERKLU);
        break;
    case SolverType::ESPQR:
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCQR);
        PCFactorSetMatSolverType(pc, MATSOLVERSPQR);
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
        KSPSetType(ksp,KSPIBCGS);
        KSPGetPC(ksp, &pc);
        PCSetType(pc,PCJACOBI);
        KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,200);
        break;

    default:
        PanicButton();
        break;
    }

    // KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    auto * pmat = dynamic_cast<PETScMatrix*> (fAnalysis->GlobalMatrix());

    KSPSolve(ksp,pmat->Rhs(),pmat->Solution());
    // KSPGetTotalIterations(ksp, &iterations); 

    // VecView(fAnalysis->Solution(),PETSC_VIEWER_STDOUT_WORLD);

}

#endif