#ifndef TENSOR_H
#define TENSOR_H

#include "DenseEigen.h"
#include "Tensor4D.h"

#define XX 0
#define YY 1
#define ZZ 2
#define XY 5
#define XZ 4
#define YZ 3

class Tensor4D;

class Tensor{
public:
    VecDouble    fData;

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
    void Identity();
    Tensor Normalized();
    Tensor NormalizedDeviatory();
    Tensor Multiply(MatrixDouble &mat);

    Tensor operator+(const Tensor &sum) const;
    const Tensor & operator+=(const Tensor &sum);

    Tensor operator-(const Tensor &sum) const;
    const Tensor & operator-=(const Tensor &sum);

    Tensor operator*(const double &multipl) const;
    const Tensor & operator*=(const double &multipl);

    Tensor operator/(const double &multipl) const;
    const Tensor & operator/=(const double &multipl);    
    
    void Zero();

    VecDouble Eigenvalues();

    /// @brief Computes Eigenvalues and eigenvectors
    /// @return 
    void SpectralDecomposition(VecDouble &eigenvalues, std::vector<MatrixDouble> &eigenprojections);


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

    MatrixDouble MatrixForm();
    
    MatrixDouble TensorProduct(Tensor &tensor);

};





#endif