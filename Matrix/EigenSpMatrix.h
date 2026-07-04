#ifndef EIGEN_SPMATRIX_H
#define EIGEN_SPMATRIX_H
#include "MatrixType.h"
#include "SparseEigen.h"
#include "DataTypes.h"

class EigenSpMatrix : public MatrixType
{
protected:
    SparseMat     fMatrix;
    MatrixDouble  fSolution;
    MatrixDouble  fRhs;

public:
    EigenSpMatrix(int rows, int cols);

    ~EigenSpMatrix();

    void AddValueMatrix(int &row, int &col, REAL &val);
    void PutValueMatrix(int &row, int &col, REAL &val);
    REAL GetValueMatrix(int &row, int &col);

    void AddValueRhs(int &row, REAL &val);
    void PutValueRhs(int &row, REAL &val);
    REAL GetValueRhs(int &row);

    void AddValueSolution(int &row, REAL &val);
    void PutValueSolution(int &row, REAL &val);
    REAL GetValueSolution(int &row);
    
    virtual SparseMat &Matrix(){return fMatrix;}
    virtual MatrixDouble &Rhs(){return fRhs;}
    virtual MatrixDouble &Solution(){return fSolution;}
    
    void ZeroMatrix() override;
    void ZeroRhs() override;
    void ZeroSolution() override;

    void PrintMatrix() override;
    void PrintRhs() override;
    void PrintSolution() override;

    REAL SolutionNorm() override;

};




#endif
