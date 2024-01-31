#include "Tensor.h"
#include "PanicButton.h"

Tensor::Tensor(){
    std::cout << "Please provide a tensor" << std::endl;
    PanicButton();
}

Tensor::Tensor(MatrixDouble &tensor){
    fTensor = tensor;
    fDimension = fTensor.rows();
#ifdef DEBUG_BUILD
    if (fTensor.rows() != fTensor.cols()){
        std::cout << "Please provide a symmetric tensor" << std::endl;
        PanicButton();
    }
#endif
}

void Tensor::ComputeHydrostatic(){
    fHydrostatic.resize(fDimension,fDimension);
    fHydrostatic.setIdentity();
    fHydrostatic *= I1() / 3.;
}


MatrixDouble Tensor::Hydrostatic(){
    if (fHydrostatic.rows() != 0){
        return fHydrostatic;
    } else {
        ComputeHydrostatic();
    }
    return fHydrostatic;
}

void Tensor::ComputeDeviatory(){
    fDeviatory = fTensor - Hydrostatic();
}

MatrixDouble Tensor::Deviatory(){
    if (fDeviatory.rows() != 0){
        return fDeviatory;
    } else {
        ComputeDeviatory();
    }
    return fDeviatory;
}

double Tensor::I1(){
    return fTensor.trace();
}

double Tensor::I2(){
    double i2= 0.;
    for (int i = 0; i < fDimension; i++){
        for (int j = 0; j < fDimension; j++){
            i2 += 0.5 * (fTensor(i,i) * fTensor(j,j) - fTensor(i,j) * fTensor(j,i));
        }
    }
    return i2;
}

double Tensor::I3(){
    return fTensor.determinant();
}

double Tensor::J1(){
    Deviatory();
    double j1 = fDeviatory.trace();
#ifdef DEBUG_BUILD
    if (fabs(j1) > 1.e-10) {
        PanicButton();
    }
#endif
    return j1;
}

double Tensor::J2(){
    double i1 = I1();
    double j2 = i1*i1 / 3. - I2();
    return j2;
}


double Tensor::J3(){
    Deviatory();
    return fDeviatory.determinant();
}