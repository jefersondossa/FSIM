#ifndef PETSC_MATRIX_H
#define PETSC_MATRIX_H
#include "MatrixType.h"
#include "PETScDataTypes.h"

#ifdef HAS_PETSC
enum PETScMatType {ESeq, ESeqSym, EAij};
#include <petscerror.h>

class PETScMatrix : public MatrixType
{
    
private:
    PETScMat fMatrix;
    PETScVec fRhs;
    PETScVec fSolution;
    PETScVec fSolutionAll;

public:
    
    PETScMatrix(int rows, int cols, PETScMatType mtype = ESeq);
    ~PETScMatrix();

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

    PETScMat &Matrix() {
        return fMatrix;
    }
    PETScVec &Rhs() {
        return fRhs;
    }
    PETScVec &Solution() {
        return fSolution;
    }

    void ExpandSolution() override;
    void ClearSolution() override;
    REAL SolutionNorm() override;
};

#endif

#endif