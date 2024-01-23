#ifndef EIGEN_SPMATRIX_H
#define EIGEN_SPMATRIX_H
#include "MatrixType.h"
#include "SparseEigen.h"

class EigenSpMatrix : public MatrixType
{
protected:
    SparseMat *fMatrix;

public:
    EigenSpMatrix(int64_t rows, int64_t cols);

    ~EigenSpMatrix();

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
