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
    virtual ~MatrixType();

    virtual void AddValueMatrix(int64_t &row, int64_t &col, double &val) = 0;
    virtual void PutValueMatrix(int64_t &row, int64_t &col, double &val) = 0;
    virtual double GetValueMatrix(int64_t &row, int64_t &col) = 0;

    virtual void AddValueRhs(int64_t &row, double &val) = 0;
    virtual void PutValueRhs(int64_t &row, double &val) = 0;
    virtual double GetValueRhs(int64_t &row) = 0;

    virtual void AddValueSolution(int64_t &row, double &val) = 0;
    virtual void PutValueSolution(int64_t &row, double &val) = 0;
    virtual double GetValueSolution(int64_t &row) = 0;

    virtual void PrintMatrix(){};
    virtual void PrintRhs(){};
    virtual void PrintSolution(){};

    virtual void ZeroMatrix(){};
    virtual void ZeroRhs(){};
    virtual void ZeroSolution(){};

    virtual void MatVecAssemble(){};

    virtual void ExpandSolution(){};
    virtual void ClearSolution(){};

    virtual double SolutionNorm() = 0; 

};




#endif