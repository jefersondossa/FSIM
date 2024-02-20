#include "Tensor.h"
#include "PanicButton.h"

Tensor::Tensor(){
    fData.resize(6);
    fData.setZero();
}

Tensor::Tensor(const Tensor &tensor){
    fData = tensor.fData;
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
    fData[XY] = tensor(0,1);
    fData[YY] = tensor(1,1);
    fData[XY] = tensor(0,1);
    fData[YY] = tensor(1,1);
    fData[XZ] = tensor(0,2);
    fData[YZ] = tensor(1,2);
    fData[YZ] = tensor(2,2);
}

void Tensor::SetData(MatrixDouble &tensor){
#ifdef DEBUG_BUILD
    if (tensor.rows() != tensor.cols() || tensor.rows() == 0){
        std::cout << "Please provide a symmetric tensor" << std::endl;
        PanicButton();
    }
#endif
    fData[XX] = tensor(0,0);
    fData[XY] = tensor(0,1);
    fData[YY] = tensor(1,1);
    fData[XY] = tensor(0,1);
    fData[YY] = tensor(1,1);
    fData[XZ] = tensor(0,2);
    fData[YZ] = tensor(1,2);
    fData[YZ] = tensor(2,2);
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

const Tensor & Tensor::operator/=(const double &multipl) {
    int i;
    for (i = 0; i < 6; i++)fData[i] /= multipl;
    return *this;
}

Tensor Tensor::operator/(const double &multipl) const {
    Tensor temp(*this);
    return temp /= multipl;
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

const Tensor & Tensor::operator-=(const Tensor &sum) {
    int i;
    for (i = 0; i < 6; i++)fData[i] -= sum.fData[i];
    return *this;
}

Tensor Tensor::operator-(const Tensor &sum) const {
    Tensor temp(*this);
    return temp -= sum;
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
    VecDouble eigval(3);
    eigval.setZero();
    
    MatrixDouble aux(3,3);
    aux(0,0) = fData[XX];
    aux(1,1) = fData[YY];
    aux(2,2) = fData[ZZ];
    aux(0,1) = aux(1,0) = fData[XY];
    aux(0,2) = aux(2,0) = fData[XZ];
    aux(1,2) = aux(2,1) = fData[YZ];

    EigenSolver<MatrixDouble> solver(aux,EigenvaluesOnly);
    eigval = solver.eigenvalues().real();
    
    return eigval;
}

Tensor Tensor::Multiply(MatrixDouble &mat){
    Tensor temp(*this);
#ifdef DEBUG_BUILD
    if(mat.rows() != 6 || mat.cols() != 6){
        std::cout << "Matrix with wrong size\n";
        PanicButton();
    }
#endif

    VecDouble aux = mat * fData;
    temp.fXX() = aux[XX];
    temp.fYY() = aux[YY];
    temp.fZZ() = aux[ZZ];
    temp.fXY() = aux[XY];
    temp.fXZ() = aux[XZ];
    temp.fYZ() = aux[YZ];
    
    return temp;
}

MatrixDouble Tensor::MatrixForm(){
    MatrixDouble mat(3,3);
    mat.setZero();

    mat(0,0) = fData[XX];
    mat(1,1) = fData[YY];
    mat(2,2) = fData[ZZ];
    mat(0,1) = mat(1,0) = fData[XY];
    mat(0,2) = mat(2,0) = fData[XZ];
    mat(1,2) = mat(2,1) = fData[YZ];

    return mat;
}

MatrixDouble Tensor::TensorProduct(Tensor &tensor){
    MatrixDouble tensor1 = this->MatrixForm();
    MatrixDouble tensor2 = tensor.MatrixForm();
    //Based on this reference: https://wiki.seg.org/wiki/Voigt_notation
    
    MatrixDouble result(6,6);
    result.setZero();
    //Diagonal
    result(0,0) = tensor1(0,0) * tensor2(0,0);
    result(1,1) = tensor1(1,1) * tensor2(1,1);
    result(2,2) = tensor1(2,2) * tensor2(2,2);
    result(3,3) = tensor1(1,2) * tensor2(1,2);
    result(4,4) = tensor1(0,2) * tensor2(0,2);
    result(5,5) = tensor1(0,1) * tensor2(0,1);
    //Off diagonal
    //1st line
    result(0,1) = result(1,0) = tensor1(0,0) * tensor2(1,1);
    result(0,2) = result(2,0) = tensor1(0,0) * tensor2(2,2);
    result(0,3) = result(3,0) = tensor1(0,0) * tensor2(1,2);
    result(0,4) = result(4,0) = tensor1(0,0) * tensor2(0,2);
    result(0,5) = result(5,0) = tensor1(0,0) * tensor2(0,1);
    //2nd line
    result(1,2) = result(2,1) = tensor1(1,1) * tensor2(2,2);
    result(1,3) = result(3,1) = tensor1(1,1) * tensor2(1,2);
    result(1,4) = result(4,1) = tensor1(1,1) * tensor2(0,2);
    result(1,5) = result(5,1) = tensor1(1,1) * tensor2(0,1);
    //3rd line
    result(2,3) = result(3,2) = tensor1(2,2) * tensor2(1,2);
    result(2,4) = result(4,2) = tensor1(2,2) * tensor2(0,2);
    result(2,5) = result(5,2) = tensor1(2,2) * tensor2(0,1);
    //4th line
    result(3,4) = result(4,3) = tensor1(1,2) * tensor2(0,2);
    result(3,5) = result(5,3) = tensor1(1,2) * tensor2(0,1);
    //5th line
    result(4,5) = result(5,4) = tensor1(0,2) * tensor2(0,1);

    return result;
}