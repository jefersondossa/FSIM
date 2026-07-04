#ifndef TENSOR3D_H
#define TENSOR3D_H

#include "DenseEigen.h"
#include "Tensor4D.h"

#define XX 0
#define YY 1
#define ZZ 2
#define XY 5
#define XZ 4
#define YZ 3

class Tensor4D;

class Tensor3D{
public:
    VecDouble    fData;

public:
    Tensor3D();

    Tensor3D(MatrixDouble &tensor);
    Tensor3D(VecDouble &tensor);
    Tensor3D(const Tensor3D &tensor);

    ~Tensor3D() = default;

    void SetData(MatrixDouble &tensor);

    //Tensor invariants
    REAL I1() const;
    REAL I2() const;
    REAL I3() const;

    REAL J1();
    REAL J2() const;
    REAL J3();

    REAL Trace() const;

    Tensor3D Deviatory();

    Tensor3D Hydrostatic();

    REAL Determinant();

    REAL Norm();
    REAL DeviatoryNorm();

    REAL DoubleContraction(Tensor3D &t);
    void Identity();
    Tensor3D Normalized();
    Tensor3D NormalizedDeviatory();
    Tensor3D Multiply(MatrixDouble &mat);

    Tensor3D operator+(const Tensor3D &sum) const;
    const Tensor3D & operator+=(const Tensor3D &sum);

    Tensor3D operator-(const Tensor3D &sum) const;
    const Tensor3D & operator-=(const Tensor3D &sum);

    Tensor3D operator*(const REAL &multipl) const;
    const Tensor3D & operator*=(const REAL &multipl);

    Tensor3D operator/(const REAL &multipl) const;
    const Tensor3D & operator/=(const REAL &multipl);

    Tensor3D & operator=(const Tensor3D &source){
        fData = source.fData;
        return *this;
    };


    friend std::ostream &operator<<(std::ostream &out, const Tensor3D &tens){
        out << tens.fData << std::endl;
        return out;
    };
    
    void Zero();

    VecDouble Eigenvalues();
    MatrixDouble Eigenvectors();

    /// @brief Computes Eigenvalues and eigenvectors
    /// @return 
    void SpectralDecomposition(VecDouble &eigenvalues, std::vector<MatrixDouble> &eigenprojections);


    REAL &fXX(){
        return fData[XX];
    }
    REAL &fYY(){
        return fData[YY];
    }
    REAL &fZZ(){
        return fData[ZZ];
    }
    REAL &fXY(){
        return fData[XY];
    }
    REAL &fXZ(){
        return fData[XZ];
    }
    REAL &fYZ(){
        return fData[YZ];
    }

    MatrixDouble MatrixForm();

    MatrixDouble TensorProduct(Tensor3D &tensor);

};





#endif