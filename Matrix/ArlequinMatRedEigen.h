#ifndef ArlequinMatRedEigen_H
#define ArlequinMatRedEigen_H
#include "MatrixType.h"
#include "SparseEigen.h"
#include "DataTypes.h"

class ArlequinMatRedEigen : public MatrixType
{
protected:
    int       fDim0,fDim1,fDimLagrange;
    SparseMat     fK0,fK1,fL0,fL1,fE,fL0T,fL1T;
    MatrixDouble  fU0,fU1,fLambda;
    MatrixDouble  fF0,fF1,fG;

public:
    ArlequinMatRedEigen(int dim0, int dim1, int dimLagrange);

    ~ArlequinMatRedEigen();

    void AddValueMatrix(int &row, int &col, REAL &val);
    void PutValueMatrix(int &row, int &col, REAL &val);
    REAL GetValueMatrix(int &row, int &col);

    void AddValueRhs(int &row, REAL &val);
    void PutValueRhs(int &row, REAL &val);
    REAL GetValueRhs(int &row);

    void AddValueSolution(int &row, REAL &val);
    void PutValueSolution(int &row, REAL &val);
    REAL GetValueSolution(int &row);
    
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

    REAL SolutionNorm() override;

};




#endif
