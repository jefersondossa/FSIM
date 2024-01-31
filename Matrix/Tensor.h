#ifndef TENSOR_H
#define TENSOR_H

#include "DenseEigen.h"

class Tensor{
protected:
    MatrixDouble fTensor;
    MatrixDouble fDeviatory;
    MatrixDouble fHydrostatic;
    double fDimension;

public:
    Tensor();

    Tensor(MatrixDouble &tensor);
    
    ~Tensor() = default;

    double I1();
    double I2();
    double I3();

    double J1();
    double J2();
    double J3();

    void ComputeDeviatory();
    void ComputeHydrostatic();

    MatrixDouble Deviatory();

    MatrixDouble Hydrostatic();


};





#endif