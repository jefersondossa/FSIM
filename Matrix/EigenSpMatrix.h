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

    void AddValueMatrix(int &row, int &col, double &val);
    void PutValueMatrix(int &row, int &col, double &val);
    double GetValueMatrix(int &row, int &col);

    void AddValueRhs(int &row, double &val);
    void PutValueRhs(int &row, double &val);
    double GetValueRhs(int &row);

    void AddValueSolution(int &row, double &val);
    void PutValueSolution(int &row, double &val);
    double GetValueSolution(int &row);
    
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
