#ifndef PETSC_MATRIX_H
#define PETSC_MATRIX_H
#include "MatrixType.h"
#include "PETScDataTypes.h"

#ifdef HAS_PETSC

class PETScMatrix : public MatrixType
{
private:
    PETScMat *fMatrix;
    PETScVec *fRhs;
    PETScVec *fSolution;

public:
    PETScMatrix(int64_t rows, int64_t cols);
    ~PETScMatrix();

    void AddValue(int64_t row, int64_t col, double val) override{
        PanicButton();
    };
    void PutValue(int64_t row, int64_t col, double val) override{
        PanicButton();
    };
    double GetValue(int64_t row, int64_t col) override{
        PanicButton();
    };
};

#endif

#endif