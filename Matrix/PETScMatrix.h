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
    
    PETScMatrix(int64_t rows, int64_t cols, PETScMatType mtype = ESeq);
    ~PETScMatrix();

    void AddValueMatrix(int64_t &row, int64_t &col, double &val) override;
    void PutValueMatrix(int64_t &row, int64_t &col, double &val) override;
    double GetValueMatrix(int64_t &row, int64_t &col) override;

    void AddValueRhs(int64_t &row, double &val) override;
    void PutValueRhs(int64_t &row, double &val) override;
    double GetValueRhs(int64_t &row) override;

    void AddValueSolution(int64_t &row, double &val) override;
    void PutValueSolution(int64_t &row, double &val) override;
    double GetValueSolution(int64_t &row) override;

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
    double SolutionNorm() override;
};

#endif

#endif