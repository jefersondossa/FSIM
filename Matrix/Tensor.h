#ifndef TENSOR_H
#define TENSOR_H

#include "DenseEigen.h"

#define XX 0
#define YY 1
#define ZZ 2
#define XY 3
#define XZ 4
#define YZ 5

class Tensor{
public:
    VecDouble    fData;
    int          fDimension;
public:
    Tensor();

    Tensor(MatrixDouble &tensor);
    Tensor(const Tensor &tensor);
    
    ~Tensor() = default;

    void SetData(MatrixDouble &tensor);

    //Tensor invariants
    double I1() const;
    double I2() const;
    double I3() const;

    double J1();
    double J2() const;
    double J3();

    double Trace() const;

    Tensor Deviatory();

    Tensor Hydrostatic();

    double Determinant();

    double Norm();
    double DeviatoryNorm();

    double DoubleContraction(Tensor &t);
    Tensor TensorProduct(Tensor &t);
    void Identity();
    Tensor Normalized();
    Tensor NormalizedDeviatory();

    const Tensor & operator*=(const double &multipl);
    Tensor operator*(const double &multipl) const;
    Tensor operator+(const Tensor &sum) const;
    const Tensor & operator+=(const Tensor &sum);
    
    void Zero();

    VecDouble Eigenvalues();


    double &fXX(){
        return fData[XX];
    }
    double &fYY(){
        return fData[YY];
    }
    double &fZZ(){
        return fData[ZZ];
    }
    double &fXY(){
        return fData[XY];
    }
    double &fXZ(){
        return fData[XZ];
    }
    double &fYZ(){
        return fData[YZ];
    }
    
};





#endif