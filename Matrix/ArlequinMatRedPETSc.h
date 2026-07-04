#ifndef ArlequinMatRedPETSc_H
#define ArlequinMatRedPETSc_H
#include "MatrixType.h"
#include "PETScDataTypes.h"

#ifdef HAS_PETSC
#include <petscerror.h>
#include "PETScMatrix.h"

class ArlequinMatRedPETSc : public MatrixType
{
    
private:
    int  fDim0,fDim1,fDimLagrange;
    PETScMat fK0,fK1,fL0,fL1,fE,fL0T,fL1T;
    PETScVec fF0,fF1,fG;
    PETScVec fU0,fU1,fLambda;
    PETScVec fSolutionAll;

public:
    
    ArlequinMatRedPETSc(int dim0, int dim1, int dimLagrange, PETScMatType mtype = ESeq);
    ~ArlequinMatRedPETSc();

    void AddValueMatrix(int &row, int &col, REAL &val) override;
    void PutValueMatrix(int &row, int &col, REAL &val) override;
    REAL GetValueMatrix(int &row, int &col) override;

    void AddValueRhs(int &row, REAL &val) override;
    void PutValueRhs(int &row, REAL &val) override;
    REAL GetValueRhs(int &row) override;

    void AddValueSolution(int &row, REAL &val) override;
    void PutValueSolution(int &row, REAL &val) override;
    REAL GetValueSolution(int &row) override;

    void ZeroMatrix() override;
    void ZeroRhs() override;
    void ZeroSolution() override;

    void PrintMatrix() override;
    void PrintRhs() override;
    void PrintSolution() override;

    void MatVecAssemble() override;
    void MatAssemble() override;
    void VecAssemble() override;

    PETScMat &K0() {return fK0;}
    PETScMat &K1() {return fK1;}
    PETScMat &L0() {return fL0;}
    PETScMat &L1() {return fL1;}
    PETScMat &L0T() {return fL0T;}
    PETScMat &L1T() {return fL1T;}
    PETScMat &E() {return fE;}

    PETScVec &F0() {return fF0;}
    PETScVec &F1() {return fF1;}
    PETScVec &G() {return fG;}

    PETScVec &U0() {return fU0;}
    PETScVec &U1() {return fU1;}
    PETScVec &Lambda() {return fLambda;}


    void ExpandSolution() override;
    void ClearSolution() override;
    REAL SolutionNorm() override;
};

#endif

#endif