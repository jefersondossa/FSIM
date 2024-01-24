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
    EigenSpMatrix(int64_t rows, int64_t cols);

    ~EigenSpMatrix();

    void AddValueMatrix(int64_t &row, int64_t &col, double &val);
    void PutValueMatrix(int64_t &row, int64_t &col, double &val);
    double GetValueMatrix(int64_t &row, int64_t &col);

    void AddValueRhs(int64_t &row, double &val);
    void PutValueRhs(int64_t &row, double &val);
    double GetValueRhs(int64_t &row);

    void AddValueSolution(int64_t &row, double &val);
    void PutValueSolution(int64_t &row, double &val);
    double GetValueSolution(int64_t &row);
    
    virtual SparseMat &Matrix(){return fMatrix;}
    virtual MatrixDouble &Rhs(){return fRhs;}
    virtual MatrixDouble &Solution(){return fSolution;}
    
    void ZeroMatrix() override;
    void ZeroRhs() override;
    void ZeroSolution() override;

    void PrintMatrix() override;
    void PrintRhs() override;
    void PrintSolution() override;

    double SolutionNorm() override;

};




#endif
