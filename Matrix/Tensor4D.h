#ifndef TENSOR4D_H
#define TENSOR4D_H

#include "DenseEigen.h"
#include "Tensor3D.h"

class Tensor4D
{
private:
    MatrixDouble fData;
public:
    Tensor4D();
    ~Tensor4D();
    
};




#endif
