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

    void AddValueMatrix(int64_t &row, int64_t &col, double &val) override{
        PanicButton();
    };
    void PutValueMatrix(int64_t &row, int64_t &col, double &val) override{
        PanicButton();
    };
    double GetValueMatrix(int64_t &row, int64_t &col) override{
        PanicButton();
    };

    void AddValueRhs(int64_t &row, double &val) override{
        PanicButton();
    };
    void PutValueRhs(int64_t &row, double &val) override{
        PanicButton();
    };
    double GetValueRhs(int64_t &row) override{
        PanicButton();
    };

    void AddValueSolution(int64_t &row, double &val) override{
        PanicButton();
    };
    void PutValueSolution(int64_t &row, double &val) override{
        PanicButton();
    };
    double GetValueSolution(int64_t &row) override{
        PanicButton();
    };
    
    virtual SparseMat* Matrix(){return fMatrix;}
    double SolutionNorm() override{return 0;}

};




#endif
