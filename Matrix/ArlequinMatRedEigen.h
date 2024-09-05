#ifndef ArlequinMatRedEigen_H
#define ArlequinMatRedEigen_H
#include "MatrixType.h"
#include "SparseEigen.h"
#include "DataTypes.h"

class ArlequinMatRedEigen : public MatrixType
{
protected:
    int64_t       fDim0,fDim1,fDimLagrange;
    SparseMat     fK0,fK1,fL0,fL1,fE,fL0T,fL1T;
    MatrixDouble  fU0,fU1,fLambda;
    MatrixDouble  fF0,fF1,fG;

public:
    ArlequinMatRedEigen(int64_t dim0, int64_t dim1, int64_t dimLagrange);

    ~ArlequinMatRedEigen();

    void AddValueMatrix(int64_t &row, int64_t &col, double &val);
    void PutValueMatrix(int64_t &row, int64_t &col, double &val);
    double GetValueMatrix(int64_t &row, int64_t &col);

    void AddValueRhs(int64_t &row, double &val);
    void PutValueRhs(int64_t &row, double &val);
    double GetValueRhs(int64_t &row);

    void AddValueSolution(int64_t &row, double &val);
    void PutValueSolution(int64_t &row, double &val);
    double GetValueSolution(int64_t &row);
    
    virtual SparseMat &K0(){return fK0;}
    virtual SparseMat &K1(){return fK1;}
    virtual SparseMat &L0(){return fL0;}
    virtual SparseMat &L1(){return fL1;}
    virtual SparseMat &L0T(){return fL0T;}
    virtual SparseMat &L1T(){return fL1T;}
    virtual SparseMat &E(){return fE;}
    virtual MatrixDouble &F0(){return fF0;}
    virtual MatrixDouble &F1(){return fF1;}
    virtual MatrixDouble &G(){return fG;}
    virtual MatrixDouble &U0(){return fU0;}
    virtual MatrixDouble &U1(){return fU1;}
    virtual MatrixDouble &Lambda(){return fLambda;}
    
    void ZeroMatrix() override;
    void ZeroRhs() override;
    void ZeroSolution() override;

    void PrintMatrix() override;
    void PrintRhs() override;
    void PrintSolution() override;

    double SolutionNorm() override;

};




#endif
