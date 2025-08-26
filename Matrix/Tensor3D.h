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
    double I1() const;
    double I2() const;
    double I3() const;

    double J1();
    double J2() const;
    double J3();

    double Trace() const;

    Tensor3D Deviatory();

    Tensor3D Hydrostatic();

    double Determinant();

    double Norm();
    double DeviatoryNorm();

    double DoubleContraction(Tensor3D &t);
    void Identity();
    Tensor3D Normalized();
    Tensor3D NormalizedDeviatory();
    Tensor3D Multiply(MatrixDouble &mat);

    Tensor3D operator+(const Tensor3D &sum) const;
    const Tensor3D & operator+=(const Tensor3D &sum);

    Tensor3D operator-(const Tensor3D &sum) const;
    const Tensor3D & operator-=(const Tensor3D &sum);

    Tensor3D operator*(const double &multipl) const;
    const Tensor3D & operator*=(const double &multipl);

    Tensor3D operator/(const double &multipl) const;
    const Tensor3D & operator/=(const double &multipl);

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

    MatrixDouble TensorProduct(Tensor3D &tensor);

};





#endif