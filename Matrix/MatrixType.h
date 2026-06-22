#ifndef MATRIX_TYPE
#define MATRIX_TYPE
#include "SparseEigen.h"
#include "PETScDataTypes.h"
#include "PanicButton.h"

class MatrixType
{
protected:
    int fRows;
    int fCols;
    
public:
    MatrixType(int rows, int cols);
    virtual ~MatrixType();

    virtual void AddValueMatrix(int &row, int &col, double &val) = 0;
    virtual void PutValueMatrix(int &row, int &col, double &val) = 0;
    virtual double GetValueMatrix(int &row, int &col) = 0;

    virtual void AddValueRhs(int &row, double &val) = 0;
    virtual void PutValueRhs(int &row, double &val) = 0;
    virtual double GetValueRhs(int &row) = 0;

    virtual void AddValueSolution(int &row, double &val) = 0;
    virtual void PutValueSolution(int &row, double &val) = 0;
    virtual double GetValueSolution(int &row) = 0;

    virtual void PrintMatrix() = 0;
    virtual void PrintRhs() = 0;
    virtual void PrintSolution() = 0;

    virtual void ZeroMatrix() = 0;
    virtual void ZeroRhs() = 0;
    virtual void ZeroSolution() = 0;

    virtual void MatVecAssemble(){};
    virtual void MatAssemble(){};
    virtual void VecAssemble(){};

    virtual void ExpandSolution(){};
    virtual void ClearSolution(){};

    virtual double SolutionNorm() = 0; 

};




#endif