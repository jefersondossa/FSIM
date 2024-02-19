#include "Tensor.h"
#include "PanicButton.h"

Tensor::Tensor(){
    fData.resize(6);
    fData.setZero();
}

Tensor::Tensor(const Tensor &tensor){
    fData = tensor.fData;
    fDimension = tensor.fDimension;
}

Tensor::Tensor(MatrixDouble &tensor){
    fData.resize(6);
    fData.setZero();
#ifdef DEBUG_BUILD
    if (tensor.rows() != tensor.cols() || tensor.rows() == 0){
        std::cout << "Please provide a symmetric tensor" << std::endl;
        PanicButton();
    }
#endif
    fData[XX] = tensor(0,0);
    fDimension = 1;
    if (tensor.rows() == 2){
        fDimension = 2;
        fData[XY] = tensor(0,1);
        fData[YY] = tensor(1,1);
    } else if (tensor.rows() == 3){
        fDimension = 3;
        fData[XY] = tensor(0,1);
        fData[YY] = tensor(1,1);
        fData[XZ] = tensor(0,2);
        fData[YZ] = tensor(1,2);
        fData[YZ] = tensor(2,2);
    }
}

void Tensor::SetData(MatrixDouble &tensor){
#ifdef DEBUG_BUILD
    if (tensor.rows() != tensor.cols() || tensor.rows() == 0){
        std::cout << "Please provide a symmetric tensor" << std::endl;
        PanicButton();
    }
#endif
    fData[XX] = tensor(0,0);
    fDimension = 1;
    if (tensor.rows() == 2){
        fDimension = 2;
        fData[XY] = tensor(0,1);
        fData[YY] = tensor(1,1);
    } else if (tensor.rows() == 3){
        fDimension = 3;
        fData[XY] = tensor(0,1);
        fData[YY] = tensor(1,1);
        fData[XZ] = tensor(0,2);
        fData[YZ] = tensor(1,2);
        fData[YZ] = tensor(2,2);
    }
}

void Tensor::Zero(){
    fData.setZero();
}

Tensor Tensor::Hydrostatic(){
    Tensor fHydrostatic(*this);
    fHydrostatic.Zero();
    double val = I1() / 3.;
    fHydrostatic.fData[XX] = val;
    fHydrostatic.fData[YY] = val;
    fHydrostatic.fData[ZZ] = val;
    return fHydrostatic;
}


Tensor Tensor::Deviatory(){ 
    Tensor fDeviatory(*this);
    double val = I1() / 3.;
    fDeviatory.fData[XX] -= val;
    fDeviatory.fData[YY] -= val;
    fDeviatory.fData[ZZ] -= val;
    return fDeviatory;
}

double Tensor::I1() const{
    return fData[XX] + fData[YY] + fData[ZZ];
}

double Tensor::I2() const{
    return -(fData[XY] * fData[XY] +
             fData[XZ] * fData[XZ] +
             fData[YZ] * fData[YZ])
           +(fData[XX] * fData[YY] +
             fData[YY] * fData[ZZ] +
             fData[XX] * fData[ZZ]);
}

double Tensor::I3() const{
    return fData[XX] * fData[YY] * fData[ZZ]
         +(fData[XY] * fData[XZ] * fData[YZ]) * 2.
         -(fData[XX] * fData[YZ] * fData[YZ] +
           fData[YY] * fData[XZ] * fData[XZ] +
           fData[ZZ] * fData[XY] * fData[XY]);
}

double Tensor::J1(){
    auto dev = Deviatory();
    return dev.fData[XX] + dev.fData[YY] + dev.fData[ZZ];
}

double Tensor::J2() const{
    double i1 = I1();
    double j2 = i1*i1 / 3. - I2();
    return j2;
}

double Tensor::Determinant(){
    return fData[XX] * fData[YY] * fData[ZZ] + fData[XY] * fData[XZ] * fData[YZ]*2. - fData[XZ] * fData[YY] * fData[XZ] -
           fData[XY] * fData[XY] * fData[ZZ] - fData[YZ] * fData[YZ] * fData[XX];
}


double Tensor::J3(){
    auto fDeviatory = Deviatory();
    double det = fDeviatory.fData[XX] * fDeviatory.fData[YY] * fDeviatory.fData[ZZ] +
                 fDeviatory.fData[XY] * fDeviatory.fData[XZ] * fDeviatory.fData[YZ]*2. - 
                 fDeviatory.fData[XZ] * fDeviatory.fData[YY] * fDeviatory.fData[XZ] -
                 fDeviatory.fData[XY] * fDeviatory.fData[XY] * fDeviatory.fData[ZZ] - 
                 fDeviatory.fData[YZ] * fDeviatory.fData[YZ] * fDeviatory.fData[XX];
    return det;
}

double Tensor::Norm(){
    double norm = 0.;
    for (unsigned int i = 0; i < 6; i++) {
        norm += fData[i] * fData[i];
    }
    norm += fData[XY] * fData[XY];
    norm += fData[XZ] * fData[XZ];
    norm += fData[YZ] * fData[YZ];
    return sqrt(norm);
}

double Tensor::DeviatoryNorm(){
    auto fDeviatory = Deviatory();
    double norm = 0.;
    for (unsigned int i = 0; i < 6; i++) {
        norm += fDeviatory.fData[i] * fDeviatory.fData[i];
    }
    norm += fDeviatory.fData[XY] * fDeviatory.fData[XY];
    norm += fDeviatory.fData[XZ] * fDeviatory.fData[XZ];
    norm += fDeviatory.fData[YZ] * fDeviatory.fData[YZ];
    return sqrt(norm);
}

double Tensor::DoubleContraction(Tensor &t){
    return fData[XX]*t.fData[XX] + fData[YY]*t.fData[YY] + fData[ZZ]*t.fData[ZZ] + 
           2.*(fData[XY]*t.fData[XY] + fData[XZ]*t.fData[XZ] + fData[YZ]*t.fData[YZ]);
}

Tensor Tensor::TensorProduct(Tensor&t){
    PanicButton();
}

void Tensor::Identity(){
    fData.setZero();
    fData[XX] = 1.;
    fData[YY] = 1.;
    fData[ZZ] = 1.;
}

const Tensor & Tensor::operator*=(const double &multipl) {
    int i;
    for (i = 0; i < 6; i++)fData[i] *= multipl;
    return *this;
}

Tensor Tensor::operator*(const double &multipl) const {
    Tensor temp(*this);
    return temp *= multipl;
}

const Tensor & Tensor::operator+=(const Tensor &sum) {
    int i;
    for (i = 0; i < 6; i++)fData[i] += sum.fData[i];
    return *this;
}

Tensor Tensor::operator+(const Tensor &sum) const {
    Tensor temp(*this);
    return temp += sum;
}

Tensor Tensor::Normalized(){
    Tensor temp(*this);
    temp *= 1./Norm();
    return temp;
}

Tensor Tensor::NormalizedDeviatory(){
    Tensor temp = Deviatory();
    temp *= 1./DeviatoryNorm();
    return temp;
}

double Tensor::Trace() const{
    return fData[XX] + fData[YY] + fData[ZZ];
}

VecDouble Tensor::Eigenvalues(){
    VecDouble eigval(fDimension);
    eigval.setZero();
    switch (fDimension)
    {
    case 1:
        eigval[0] = fData[XX];
        break;
    case 2:
        { 
            //Using Mohr's circle relations
            double radius = sqrt((0.25*(fData[XX]-fData[YY])*(fData[XX]-fData[YY]))+fData[XY]*fData[XY]);
            double sigmed = 0.5*(fData[XX]+fData[YY]);
            eigval[0] = sigmed + radius;
            eigval[1] = sigmed - radius;
        }
        break;
    case 3:
        {
            MatrixDouble aux(3,3);
            aux(0,0) = fData[XX];
            aux(1,1) = fData[YY];
            aux(2,2) = fData[ZZ];
            aux(0,1) = aux(1,0) = fData[XY];
            aux(0,2) = aux(2,0) = fData[XZ];
            aux(1,2) = aux(2,1) = fData[YZ];

            EigenSolver<MatrixDouble> solver(aux,EigenvaluesOnly);
            eigval = solver.eigenvalues().real();
        }
    
    default:
        PanicButton();
        break;
    }

    return eigval;
}