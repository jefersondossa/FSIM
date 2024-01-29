#include "PETScMatrix.h"

#ifdef HAS_PETSC

PETScMatrix::PETScMatrix(int64_t rows, int64_t cols, PETScMatType mtype) : MatrixType(rows,cols){
    switch (mtype)
    {
    case ESeq:
        MatCreateSeqAIJ(PETSC_COMM_WORLD, rows, cols, 100,NULL,&fMatrix);
        break;
    case ESeqSym:
        MatCreateSeqAIJ(PETSC_COMM_WORLD, rows, cols, 100,NULL,&fMatrix);
        MatSetOption(fMatrix, MAT_SYMMETRIC, PETSC_TRUE);
        break;
    case EAij:
        MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                        rows, cols,100,NULL,300,NULL,&fMatrix); 
        break;
    
    default:
        PanicButton();
        break;
    }
    //Create PETSc vectors
    VecCreate(PETSC_COMM_WORLD,&fRhs);
    VecSetSizes(fRhs,PETSC_DECIDE,rows);
    
    VecSetFromOptions(fRhs);
    VecDuplicate(fRhs,&fSolution);
}

PETScMatrix::~PETScMatrix(){
    MatDestroy(&fMatrix);
    VecDestroy(&fRhs);
    VecDestroy(&fSolution);
}


void PETScMatrix::AddValueMatrix(int64_t &row, int64_t &col, double &val) {
    MatSetValues(fMatrix,1,&row,1,&col,&val,ADD_VALUES);
};

void PETScMatrix::PutValueMatrix(int64_t &row, int64_t &col, double &val) {
    MatSetValues(fMatrix,1,&row,1,&col,&val,INSERT_VALUES);
};

double PETScMatrix::GetValueMatrix(int64_t &row, int64_t &col) {
    std::cout << "It may need a scatter context. Please check PETSc manual and implement it. \n";
    PanicButton();
};

void PETScMatrix::ZeroMatrix(){
    MatZeroEntries(fMatrix);
}

void PETScMatrix::ZeroRhs(){
    VecZeroEntries(fRhs);
}

void PETScMatrix::ZeroSolution(){
    VecZeroEntries(fSolution);
}

void PETScMatrix::PrintMatrix(){
    MatView(fMatrix,PETSC_VIEWER_STDOUT_WORLD);
    // MatView(fAnalysis->Stiffness(),PETSC_VIEWER_DRAW_WORLD);
}

void PETScMatrix::PrintRhs(){
    VecView(fRhs,PETSC_VIEWER_STDOUT_WORLD);
}

void PETScMatrix::PrintSolution(){
    VecView(fSolution,PETSC_VIEWER_STDOUT_WORLD);
}

void PETScMatrix::MatVecAssemble(){
    MatAssemblyBegin(fMatrix,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fMatrix,MAT_FINAL_ASSEMBLY);
    
    VecAssemblyBegin(fRhs);
    VecAssemblyEnd(fRhs);
}

void PETScMatrix::AddValueRhs(int64_t &row, double &val) {
    VecSetValues(fRhs, 1, &row, &val, ADD_VALUES);
};

void PETScMatrix::PutValueRhs(int64_t &row, double &val) {
    VecSetValues(fRhs, 1, &row, &val, INSERT_VALUES);
};

double PETScMatrix::GetValueRhs(int64_t &row) {
    std::cout << "This vector need a scatter context to get values. Please implement it. \n";
    PanicButton();
};

void PETScMatrix::AddValueSolution(int64_t &row, double &val) {
    VecSetValues(fSolution, 1, &row, &val, ADD_VALUES);
};

void PETScMatrix::PutValueSolution(int64_t &row, double &val) {
    VecSetValues(fSolution, 1, &row, &val, INSERT_VALUES);
};

double PETScMatrix::GetValueSolution(int64_t &row) {
    int64_t Ione = 1;
    double val;
    VecGetValues(fSolutionAll, Ione, &row, &val);
    return val;
};

void PETScMatrix::ExpandSolution(){
    VecScatter        ctx;

    //Gathers the solution vector to the master process
    VecScatterCreateToAll(fSolution, &ctx, &fSolutionAll);
    VecScatterBegin(ctx, fSolution, fSolutionAll, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, fSolution, fSolutionAll, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&ctx);
}

void PETScMatrix::ClearSolution(){
    VecDestroy(&fSolutionAll);
}

double PETScMatrix::SolutionNorm(){
    double val;
    VecNorm(fSolution,NORM_2,&val);
    return val;
}

#endif