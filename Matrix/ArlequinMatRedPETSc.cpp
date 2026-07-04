#include "ArlequinMatRedPETSc.h"

#ifdef HAS_PETSC

ArlequinMatRedPETSc::ArlequinMatRedPETSc(int dim0, int dim1, int dimLagrange, PETScMatType mtype) : MatrixType(dim0+dim1+dimLagrange,dim0+dim1+dimLagrange){
    
    fDim0 = dim0;
    fDim1 = dim1;
    fDimLagrange = dimLagrange;

    switch (mtype)
    {
    case ESeq:
        MatCreateSeqAIJ(PETSC_COMM_WORLD, dim0, dim0, 100,NULL,&fK0);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, dim1, dim1, 100,NULL,&fK1);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, dimLagrange, dimLagrange, 100,NULL,&fE);
        MatCreateDense(PETSC_COMM_WORLD, dim0, dimLagrange, PETSC_DECIDE, PETSC_DECIDE, NULL,&fL0);
        MatCreateDense(PETSC_COMM_WORLD, dimLagrange, dim0, PETSC_DECIDE, PETSC_DECIDE, NULL,&fL0T);
        MatCreateDense(PETSC_COMM_WORLD, dim1, dimLagrange, PETSC_DECIDE, PETSC_DECIDE, NULL,&fL1);
        MatCreateDense(PETSC_COMM_WORLD, dimLagrange, dim1, PETSC_DECIDE, PETSC_DECIDE, NULL,&fL1T);
        break;
    case ESeqSym:
        MatCreateSeqAIJ(PETSC_COMM_WORLD, dim0, dim0, 100,NULL,&fK0);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, dim1, dim1, 100,NULL,&fK1);
        MatCreateSeqAIJ(PETSC_COMM_WORLD, dimLagrange, dimLagrange, 100,NULL,&fE);
        MatCreateDense(PETSC_COMM_WORLD, dim0, dimLagrange, PETSC_DECIDE, PETSC_DECIDE, NULL,&fL0);
        MatCreateDense(PETSC_COMM_WORLD, dimLagrange, dim0, PETSC_DECIDE, PETSC_DECIDE, NULL,&fL0T);
        MatCreateDense(PETSC_COMM_WORLD, dim1, dimLagrange, PETSC_DECIDE, PETSC_DECIDE, NULL,&fL1);
        MatCreateDense(PETSC_COMM_WORLD, dimLagrange, dim1, PETSC_DECIDE, PETSC_DECIDE, NULL,&fL1T);
        MatSetOption(fK0, MAT_SYMMETRIC, PETSC_TRUE);
        MatSetOption(fK1, MAT_SYMMETRIC, PETSC_TRUE);
        MatSetOption(fE, MAT_SYMMETRIC, PETSC_TRUE);
        break;
    case EAij:
        MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, dim0, dim0,100,NULL,100,NULL,&fK0);
        MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, dim1, dim1,100,NULL,100,NULL,&fK1);
        MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, dimLagrange, dimLagrange,100,NULL,100,NULL,&fE);
        MatCreateDense(PETSC_COMM_WORLD, dim0, dimLagrange, PETSC_DECIDE, PETSC_DECIDE, NULL,&fL0);
        MatCreateDense(PETSC_COMM_WORLD, dimLagrange, dim0, PETSC_DECIDE, PETSC_DECIDE, NULL,&fL0T);
        MatCreateDense(PETSC_COMM_WORLD, dim1, dimLagrange, PETSC_DECIDE, PETSC_DECIDE, NULL,&fL1);
        MatCreateDense(PETSC_COMM_WORLD, dimLagrange, dim1, PETSC_DECIDE, PETSC_DECIDE, NULL,&fL1T);
        
        break;
    
    default:
        PanicButton();
        break;
    }
    //Create PETSc vectors
    VecCreate(PETSC_COMM_WORLD,&fF0);
    VecSetSizes(fF0,PETSC_DECIDE,dim0);
    VecCreate(PETSC_COMM_WORLD,&fF1);
    VecSetSizes(fF1,PETSC_DECIDE,dim1);
    VecCreate(PETSC_COMM_WORLD,&fG);
    VecSetSizes(fG,PETSC_DECIDE,dimLagrange);
    
    VecSetFromOptions(fF0);
    VecSetFromOptions(fF1);
    VecSetFromOptions(fG);
    VecDuplicate(fF0,&fU0);
    VecDuplicate(fF1,&fU1);
    VecDuplicate(fG,&fLambda);
}

ArlequinMatRedPETSc::~ArlequinMatRedPETSc(){
    MatDestroy(&fK0);
    MatDestroy(&fK1);
    MatDestroy(&fE);
    MatDestroy(&fL0);
    MatDestroy(&fL1);
    MatDestroy(&fL0T);
    MatDestroy(&fL1T);
    VecDestroy(&fF0);
    VecDestroy(&fF1);
    VecDestroy(&fG);
    VecDestroy(&fU0);
    VecDestroy(&fU1);
    VecDestroy(&fLambda);
}


void ArlequinMatRedPETSc::AddValueMatrix(int &row, int &col, REAL &val) {

    if (row < fDim0 && col < fDim0){// Belongs to K0
        MatSetValues(fK0,1,&row,1,&col,&val,ADD_VALUES);
        return;
    } else if (row < fDim0 && col > fDim0+fDim1-1){ // Belongs to L0
        PetscInt posrow = (row);
        PetscInt poscol = (col-fDim0-fDim1);
        MatSetValues(fL0,1,&posrow,1,&poscol,&val,ADD_VALUES);
        return;
    } else if (col < fDim0 && row > fDim0+fDim1-1){// Belongs to L0T
        PetscInt posrow = (row-fDim0-fDim1);
        PetscInt poscol = (col);
        MatSetValues(fL0T,1,&posrow,1,&poscol,&val,ADD_VALUES);
        return;
    } else if (row > fDim0-1 && row < fDim0+fDim1 && col > fDim0-1 && col < fDim0+fDim1){ // Belongs to K1
        PetscInt posrow = (row-fDim0);
        PetscInt poscol = (col-fDim1);
        MatSetValues(fK1,1,&posrow,1,&poscol,&val,ADD_VALUES);
        return;
    } else if (row > fDim0-1 && row < fDim0+fDim1 && col > fDim0+fDim1-1){ // Belongs to L1
        PetscInt posrow = (row-fDim0);
        PetscInt poscol = (col-fDim0-fDim1);
        MatSetValues(fL1,1,&posrow,1,&poscol,&val,ADD_VALUES);
        return;
    } else if (col > fDim0-1 && col < fDim0+fDim1 && row > fDim0+fDim1-1){ // Belongs to L1T
        PetscInt posrow = (row-fDim0-fDim1);
        PetscInt poscol = (col-fDim0);
        MatSetValues(fL1T,1,&posrow,1,&poscol,&val,ADD_VALUES);
        return;
    } else if (row > fDim0+fDim1-1 && col > fDim0+fDim1-1){ // Belongs to E
        PetscInt posrow = (row-fDim0-fDim1);
        PetscInt poscol = (col-fDim0-fDim1);
        MatSetValues(fE,1,&posrow,1,&poscol,&val,ADD_VALUES);
        return;
    } else {
        PanicButton();
    }    
    
};

void ArlequinMatRedPETSc::PutValueMatrix(int &row, int &col, REAL &val) {
    if (row < fDim0 && col < fDim0){// Belongs to K0
        PetscInt rowpetsc = (row);
        PetscInt colpetsc = (col);
        MatSetValues(fK0,1,&rowpetsc,1,&colpetsc,&val,INSERT_VALUES);
        return;
    } else if (row < fDim0 && col > fDim0+fDim1-1){ // Belongs to L0
        PetscInt posrow = (row);
        PetscInt poscol = (col-fDim0-fDim1);
        MatSetValues(fL0,1,&posrow,1,&poscol,&val,INSERT_VALUES);
        return;
    } else if (col < fDim0 && row > fDim0+fDim1-1){// Belongs to L0T
        PetscInt posrow = (row-fDim0-fDim1);
        PetscInt poscol = (col);
        MatSetValues(fL0T,1,&posrow,1,&poscol,&val,INSERT_VALUES);
        return;
    } else if (row > fDim0-1 && row < fDim0+fDim1 && col > fDim0-1 && col < fDim0+fDim1){ // Belongs to K1
        PetscInt posrow = (row-fDim0);
        PetscInt poscol = (col-fDim1);
        MatSetValues(fK1,1,&posrow,1,&poscol,&val,INSERT_VALUES);
        return;
    } else if (row > fDim0-1 && row < fDim0+fDim1 && col > fDim0+fDim1-1){ // Belongs to L1
        PetscInt posrow = (row-fDim0);
        PetscInt poscol = (col-fDim0-fDim1);
        MatSetValues(fL1,1,&posrow,1,&poscol,&val,INSERT_VALUES);
        return;
    } else if (col > fDim0-1 && col < fDim0+fDim1 && row > fDim0+fDim1-1){ // Belongs to L1T
        PetscInt posrow = row-fDim0-fDim1;
        PetscInt poscol = col-fDim0;
        MatSetValues(fL1T,1,&posrow,1,&poscol,&val,INSERT_VALUES);
        return;
    } else if (row > fDim0+fDim1-1 && col > fDim0+fDim1-1){ // Belongs to E
        PetscInt posrow = row-fDim0-fDim1;
        PetscInt poscol = col-fDim0-fDim1;
        MatSetValues(fE,1,&posrow,1,&poscol,&val,INSERT_VALUES);
        return;
    } else {
        PanicButton();
    }
    
};

REAL ArlequinMatRedPETSc::GetValueMatrix(int &row, int &col) {
    std::cout << "It may need a scatter context. Please check PETSc manual and implement it. \n";
    PanicButton();
};

void ArlequinMatRedPETSc::ZeroMatrix(){
    MatZeroEntries(fK0);
    MatZeroEntries(fK1);
    MatZeroEntries(fL0);
    MatZeroEntries(fL0T);
    MatZeroEntries(fL1);
    MatZeroEntries(fL1T);
}

void ArlequinMatRedPETSc::ZeroRhs(){
    VecZeroEntries(fF0);
    VecZeroEntries(fF1);
    VecZeroEntries(fG);
}

void ArlequinMatRedPETSc::ZeroSolution(){
    VecZeroEntries(fU0);
    VecZeroEntries(fU1);
    VecZeroEntries(fLambda);
}

void ArlequinMatRedPETSc::PrintMatrix(){
    std::cout << "K0 = \n";
    MatView(fK0,PETSC_VIEWER_STDOUT_WORLD);
    std::cout << "K1 = \n";
    MatView(fK1,PETSC_VIEWER_STDOUT_WORLD);
    std::cout << "L0 = \n";
    MatView(fL0,PETSC_VIEWER_STDOUT_WORLD);
    std::cout << "L1 = \n";
    MatView(fL1,PETSC_VIEWER_STDOUT_WORLD);
    std::cout << "E = \n";
    MatView(fE,PETSC_VIEWER_STDOUT_WORLD);
    // MatView(fAnalysis->Stiffness(),PETSC_VIEWER_DRAW_WORLD);
}

void ArlequinMatRedPETSc::PrintRhs(){
    std::cout << "F0 = \n";
    VecView(fF0,PETSC_VIEWER_STDOUT_WORLD);
    std::cout << "F1 = \n";
    VecView(fF1,PETSC_VIEWER_STDOUT_WORLD);
    std::cout << "G = \n";
    VecView(fG,PETSC_VIEWER_STDOUT_WORLD);
}

void ArlequinMatRedPETSc::PrintSolution(){
    std::cout << "U0 = \n";
    VecView(fU0,PETSC_VIEWER_STDOUT_WORLD);
    std::cout << "U1 = \n";
    VecView(fU1,PETSC_VIEWER_STDOUT_WORLD);
    std::cout << "Lambda = \n";
    VecView(fLambda,PETSC_VIEWER_STDOUT_WORLD);
}

void ArlequinMatRedPETSc::MatVecAssemble(){
    MatAssemblyBegin(fK0,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fK0,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fK1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fK1,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fL0,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fL0,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fL0T,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fL0T,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fL1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fL1,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fL1T,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fL1T,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fE,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fE,MAT_FINAL_ASSEMBLY);
    
    VecAssemblyBegin(fF0);
    VecAssemblyEnd(fF0);
    VecAssemblyBegin(fF1);
    VecAssemblyEnd(fF1);
    VecAssemblyBegin(fG);
    VecAssemblyEnd(fG);
}

void ArlequinMatRedPETSc::MatAssemble(){
    MatAssemblyBegin(fK0,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fK0,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fK1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fK1,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fL0,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fL0,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fL0T,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fL0T,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fL1,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fL1,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fL1T,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fL1T,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fE,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fE,MAT_FINAL_ASSEMBLY);
}

void ArlequinMatRedPETSc::VecAssemble(){
    VecAssemblyBegin(fF0);
    VecAssemblyEnd(fF0);
    VecAssemblyBegin(fF1);
    VecAssemblyEnd(fF1);
    VecAssemblyBegin(fG);
    VecAssemblyEnd(fG);
}

void ArlequinMatRedPETSc::AddValueRhs(int &row, REAL &val) {
    if (row < fDim0){// Belongs to F0
        VecSetValues(fF0, 1, &row, &val, ADD_VALUES);
        return;
    } else if (row >= fDim0 && row < fDim0+fDim1){ // Belongs to F1
        PetscInt posrow = row-fDim0;
        VecSetValues(fF1, 1, &posrow, &val, ADD_VALUES);
        return;
    } else if (row >= fDim0+fDim1){ // Belongs to G
        PetscInt posrow = row-fDim0-fDim1;
        VecSetValues(fG, 1, &posrow, &val, ADD_VALUES);
        return;
    } else {
        PanicButton();
    }   
};

void ArlequinMatRedPETSc::PutValueRhs(int &row, REAL &val) {
    if (row < fDim0){// Belongs to F0
        VecSetValues(fF0, 1, &row, &val, INSERT_VALUES);
        return;
    } else if (row >= fDim0 && row < fDim0+fDim1){ // Belongs to F1
        PetscInt posrow = row-fDim0;
        VecSetValues(fF1, 1, &posrow, &val, INSERT_VALUES);
        return;
    } else if (row >= fDim0+fDim1){ // Belongs to G
        PetscInt posrow = row-fDim0-fDim1;
        VecSetValues(fG, 1, &posrow, &val, INSERT_VALUES);
        return;
    } else {
        PanicButton();
    }   
};

REAL ArlequinMatRedPETSc::GetValueRhs(int &row) {
    std::cout << "This vector need a scatter context to get values. Please implement it. \n";
    PanicButton();
};

void ArlequinMatRedPETSc::AddValueSolution(int &row, REAL &val) {
    if (row < fDim0){// Belongs to U0
        VecSetValues(fU0, 1, &row, &val, ADD_VALUES);
        return;
    } else if (row >= fDim0 && row < fDim0+fDim1){ // Belongs to U1
        int posrow = row-fDim0;
        VecSetValues(fU1, 1, &posrow, &val, ADD_VALUES);
        return;
    } else if (row >= fDim0+fDim1){ // Belongs to Lambda
        int posrow = row-fDim0-fDim1;
        VecSetValues(fLambda, 1, &posrow, &val, ADD_VALUES);
        return;
    } else {
        PanicButton();
    }   
};

void ArlequinMatRedPETSc::PutValueSolution(int &row, REAL &val) {
    if (row < fDim0){// Belongs to U0
        VecSetValues(fU0, 1, &row, &val, INSERT_VALUES);
        return;
    } else if (row >= fDim0 && row < fDim0+fDim1){ // Belongs to U1
        int posrow = row-fDim0;
        VecSetValues(fU1, 1, &posrow, &val, INSERT_VALUES);
        return;
    } else if (row >= fDim0+fDim1){ // Belongs to Lambda
        int posrow = row-fDim0-fDim1;
        VecSetValues(fLambda, 1, &posrow, &val, INSERT_VALUES);
        return;
    } else {
        PanicButton();
    }
};

REAL ArlequinMatRedPETSc::GetValueSolution(int &row) {
    int Ione = 1;
    REAL val;
    VecGetValues(fSolutionAll, Ione, &row, &val);
    return val;
};

void ArlequinMatRedPETSc::ExpandSolution(){
    // VecScatter        ctx;

    // //Gathers the solution vector to the master process
    // VecScatterCreateToAll(fSolution, &ctx, &fSolutionAll);
    // VecScatterBegin(ctx, fSolution, fSolutionAll, INSERT_VALUES, SCATTER_FORWARD);
    // VecScatterEnd(ctx, fSolution, fSolutionAll, INSERT_VALUES, SCATTER_FORWARD);
    // VecScatterDestroy(&ctx);
    Vec* conc;
    PetscInt nx = 3;
    PetscMalloc1(nx, &conc);
    conc[0]=fU0;
    conc[1]=fU1;
    conc[2]=fLambda;
    
    VecConcatenate(nx, conc, &fSolutionAll, NULL);
}

void ArlequinMatRedPETSc::ClearSolution(){
    VecDestroy(&fSolutionAll);
}

REAL ArlequinMatRedPETSc::SolutionNorm(){
    REAL val;
    VecNorm(fSolutionAll,NORM_2,&val);
    return val;
}

#endif