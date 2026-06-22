#include "ArlequinRedSolverPETSc.h"

#ifdef HAS_PETSC

ArlequinRedSolverPETSc::ArlequinRedSolverPETSc(Analysis *an) : LinearSolver(an){
    
}

ArlequinRedSolverPETSc::~ArlequinRedSolverPETSc(){
    KSPDestroy(&kspU0);
    KSPDestroy(&kspU1);
    KSPDestroy(&kspLambda);
}

void ArlequinRedSolverPETSc::Solve(){
    //Create KSP context to solve the linear system
    
    auto * pmat = dynamic_cast<ArlequinMatRedPETSc*> (fAnalysis->GlobalMatrix());
    if (!pmat) PanicButton();

    KSPCreate(PETSC_COMM_WORLD,&kspU0);
    KSPSetOperators(kspU0,pmat->K0(),pmat->K0());

    KSPCreate(PETSC_COMM_WORLD,&kspU1);
    KSPSetOperators(kspU1,pmat->K1(),pmat->K1());

    switch (fAnalysis->SType())
    {
    case SolverType::EUmfpack:
        KSPGetPC(kspU0, &pcU0);
        PCSetType(pcU0, PCLU);
        PCFactorSetMatSolverType(pcU0, MATSOLVERUMFPACK);
        KSPGetPC(kspU1, &pcU1);
        PCSetType(pcU1, PCLU);
        PCFactorSetMatSolverType(pcU1, MATSOLVERUMFPACK);
        break;
    case SolverType::EKLU:
        KSPGetPC(kspU0, &pcU0);
        PCSetType(pcU0, PCLU);
        PCFactorSetMatSolverType(pcU0, MATSOLVERKLU);
        KSPGetPC(kspU1, &pcU1);
        PCSetType(pcU1, PCLU);
        PCFactorSetMatSolverType(pcU1, MATSOLVERKLU);
        break;
    case SolverType::ESPQR:
        KSPGetPC(kspU0, &pcU0);
        PCSetType(pcU0, PCQR);
        PCFactorSetMatSolverType(pcU0, MATSOLVERSPQR);
        KSPGetPC(kspU1, &pcU1);
        PCSetType(pcU1, PCQR);
        PCFactorSetMatSolverType(pcU1, MATSOLVERSPQR);
        break;
    case SolverType::ECholmod:
        KSPGetPC(kspU0, &pcU0);
        PCSetType(pcU0, PCCHOLESKY);
        PCFactorSetMatSolverType(pcU0, MATSOLVERCHOLMOD);
        KSPGetPC(kspU1, &pcU1);
        PCSetType(pcU1, PCCHOLESKY);
        PCFactorSetMatSolverType(pcU1, MATSOLVERCHOLMOD);
        break;
    case SolverType::EMumps:
        KSPGetPC(kspU0, &pcU0);
        PCSetType(pcU0, PCLU);
        PCFactorSetMatSolverType(pcU0, MATSOLVERMUMPS);
        KSPGetPC(kspU1, &pcU1);
        PCSetType(pcU1, PCLU);
        PCFactorSetMatSolverType(pcU1, MATSOLVERMUMPS);
        break;
    case SolverType::EIterative:
        KSPSetType(kspU0,KSPFGMRES);
        KSPGetPC(kspU0, &pcU0);
        PCSetType(pcU0, PCBJACOBI);
        KSPSetTolerances(kspU0,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,200);
        KSPSetType(kspU1,KSPFGMRES);
        KSPGetPC(kspU1, &pcU1);
        PCSetType(pcU1, PCBJACOBI);
        KSPSetTolerances(kspU1,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,200);
        break;

    default:
        PanicButton();
        break;
    }

    PETScMat Aux0, Aux1, Res0, Res1;
    MatDuplicate(pmat->L0(), MAT_SHARE_NONZERO_PATTERN, &Aux0);
    MatDuplicate(pmat->L1(), MAT_SHARE_NONZERO_PATTERN, &Aux1);
    
    //Aux0 = K0^-1*L0
    //Aux1 = K1^-1*L1
    KSPMatSolve(kspU0, pmat->L0(), Aux0);
    KSPMatSolve(kspU1, pmat->L1(), Aux1);
    MatMatMult(pmat->L0T(), Aux0, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Res0);
    MatMatMult(pmat->L1T(), Aux1, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Res1);
    MatAXPY(pmat->E(), -1., Res0, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(pmat->E(), -1., Res1, DIFFERENT_NONZERO_PATTERN);

    //Sets 1 to zero diagonals
    PetscInt rows,cols;
    MatGetSize(pmat->E(), &rows, &cols);
    std::vector<PetscInt> indices;
    for (PetscInt i = 0; i < rows; i++){
        PetscScalar val;
        MatGetValue(pmat->E(), i, i, &val);
        if (fabs(val) < 1.e-10){
            indices.push_back(i);
        }
    }
    for (auto i = 0; i < indices.size(); i++){
        PetscScalar val = 1.;
        MatSetValues(pmat->E(), 1, &indices[i], 1, &indices[i], &val, INSERT_VALUES);
    }
    MatAssemblyBegin(pmat->E(),MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(pmat->E(),MAT_FINAL_ASSEMBLY);
    

    // std::cout << is[0] << std::endl;

    //Compute Lambda problem Rhs
    PETScVec aF0, aF1, aG;
    VecDuplicate(pmat->F0(), &aF0);
    VecDuplicate(pmat->F1(), &aF1);
    VecDuplicate(pmat->G(), &aG);
    VecZeroEntries(aG);

    KSPSolve(kspU0, pmat->F0(), aF0);
    KSPSolve(kspU1, pmat->F1(), aF1);
    PetscScalar mone = -1.;
    VecScale(aF0, mone);
    VecScale(aF1, mone);

    MatMultAdd(pmat->L0T(), aF0, aG, aG);
    MatMultAdd(pmat->L1T(), aF1, pmat->G(), aG);

    //Create Lambda solver
    KSPCreate(PETSC_COMM_WORLD,&kspLambda);
    KSPSetOperators(kspLambda,pmat->E(),pmat->E());
    switch (fAnalysis->SType())
    {
    case SolverType::EUmfpack:
        KSPGetPC(kspLambda, &pcLambda);
        PCSetType(pcLambda, PCLU);
        PCFactorSetMatSolverType(pcLambda, MATSOLVERUMFPACK);
        break;
    case SolverType::EKLU:
        KSPGetPC(kspLambda, &pcLambda);
        PCSetType(pcLambda, PCLU);
        PCFactorSetMatSolverType(pcLambda, MATSOLVERKLU);
        break;
    case SolverType::ESPQR:
        KSPGetPC(kspLambda, &pcLambda);
        PCSetType(pcLambda, PCQR);
        PCFactorSetMatSolverType(pcLambda, MATSOLVERSPQR);
        break;
    case SolverType::ECholmod:
        KSPGetPC(kspLambda, &pcLambda);
        PCSetType(pcLambda, PCCHOLESKY);
        PCFactorSetMatSolverType(pcLambda, MATSOLVERCHOLMOD);
        break;
    case SolverType::EMumps:
        KSPGetPC(kspLambda, &pcLambda);
        PCSetType(pcLambda, PCLU);
        PCFactorSetMatSolverType(pcLambda, MATSOLVERMUMPS);
        break;
    case SolverType::EIterative:
        KSPSetType(kspLambda,KSPFGMRES);
        KSPGetPC(kspLambda, &pcLambda);
        PCSetType(pcLambda,PCBJACOBI);
        KSPSetTolerances(kspLambda,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,200);
        break;

    default:
        PanicButton();
        break;
    }
    
    //Solve Lambda
    KSPSolve(kspLambda,aG,pmat->Lambda());
    
    //Compute U0 and U1 Rhs
    VecCopy(pmat->Lambda(),aG);
    VecScale(aG,mone);
    VecZeroEntries(aF0);
    VecZeroEntries(aF1);
    MatMultAdd(pmat->L0(), aG, pmat->F0(), aF0);
    MatMultAdd(pmat->L1(), aG, pmat->F1(), aF1);

    //Solve U0 and U1 problems
    KSPSolve(kspU0,aF0,pmat->U0());
    KSPSolve(kspU1,aF1,pmat->U1());

    // VecView(pmat->U0(),PETSC_VIEWER_STDOUT_WORLD);
    // VecView(pmat->U1(),PETSC_VIEWER_STDOUT_WORLD);
    // MatView(pmat->E(),PETSC_VIEWER_STDOUT_WORLD);
    // KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    // KSPGetTotalIterations(ksp, &iterations); 

    MatDestroy(&Aux0);
    MatDestroy(&Aux1);
    MatDestroy(&Res0);
    MatDestroy(&Res1);
    VecDestroy(&aF0);
    VecDestroy(&aF1);
    VecDestroy(&aG);

}

#endif