#ifndef MATRIX_TYPE
#define MATRIX_TYPE
#include "SparseEigen.h"
#include "PETScDataTypes.h"
#include "PanicButton.h"

class MatrixType
{
protected:
    int64_t fRows;
    int64_t fCols;
    
public:
    MatrixType(int64_t rows, int64_t cols);
    ~MatrixType();

    virtual void AddValue(int64_t row, int64_t col, double val) = 0;
    virtual void PutValue(int64_t row, int64_t col, double val) = 0;
    virtual double GetValue(int64_t row, int64_t col) = 0;

};




#endif