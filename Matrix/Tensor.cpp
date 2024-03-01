#include "Tensor.h"
#include "PanicButton.h"


bool compare_head(const double& lhs, const double& rhs)
{
    return lhs > rhs;
}


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
    fData[XZ] = tensor(0,2);
    fData[YZ] = tensor(1,2);
    fData[ZZ] = tensor(2,2);
}

Tensor::Tensor(VecDouble &tensor){
    fData.resize(6);
    fData.setZero();
#ifdef DEBUG_BUILD
    if (tensor.size() != 6){
        std::cout << "Please provide a proper tensor" << std::endl;
        PanicButton();
    }
#endif
    fData = tensor;
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
    fData[XZ] = tensor(0,2);
    fData[YZ] = tensor(1,2);
    fData[ZZ] = tensor(2,2);
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
    Tensor fDeviatory = *this;
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
    MatrixDouble aux = MatrixForm();
    EigenSolver<MatrixDouble> solver(aux,EigenvaluesOnly);
    VecDouble eigval = solver.eigenvalues().real();

    return eigval;
}

MatrixDouble Tensor::Eigenvectors(){    
    MatrixDouble aux = MatrixForm();
    
    EigenSolver<MatrixDouble> solver(aux);
    MatrixDouble eigvec = solver.eigenvectors().real();

    return eigvec;
}

void Tensor::SpectralDecomposition(VecDouble &eigenvalues, std::vector<MatrixDouble> &eigenprojections){ 
    eigenvalues.resize(3); 
    eigenprojections.resize(3);

    //A.4 Souza Neto
    double i1 = I1();
    double i2 = I2();
    double i3 = I3();
    double R = (-2.*i1*i1*i1 + 9.*i1*i2 - 27.*i3) / 54.;
    double Q = (i1*i1 - 3.*i2)/9.;
    double sqQ = sqrt(Q);
    double valcos = R/(Q*sqQ);
    if (valcos < -1.) valcos = -1.;
    if (valcos >  1.) valcos = 1.;
    double theta = acos(valcos);
    if (fabs(i2) < 1.e-10){
        eigenvalues.setZero();
        return;
    } else {
        eigenvalues[0] = -2. * sqQ * cos(theta/3.) + i1/3.; 
        eigenvalues[1] = -2. * sqQ * cos((theta+2.*M_PI)/3.) + i1/3.; 
        eigenvalues[2] = -2. * sqQ * cos((theta-2.*M_PI)/3.) + i1/3.;
    }

    MatrixDouble aux = MatrixForm();
    

    // If all eigenvalues are different
    MatrixDouble ident(3,3);
    ident.setIdentity();
    if (fabs(eigenvalues[0]-eigenvalues[1]) > 1.e-10 && fabs(eigenvalues[0]-eigenvalues[2]) > 1.e-10 && fabs(eigenvalues[1]-eigenvalues[2]) > 1.e-10){
        for (size_t i = 0; i < 3; i++){
            eigenprojections[i] = eigenvalues[i]/(2.*eigenvalues[i]*eigenvalues[i]*eigenvalues[i]-i1*eigenvalues[i]*eigenvalues[i]+i3) *(aux*aux - (i1-eigenvalues[i])*aux + (i3/eigenvalues[i])*ident);
        }
#ifdef DEBUG_BUILD
    if (std::isnan(eigenprojections[0].norm())||std::isnan(eigenprojections[1].norm())||std::isnan(eigenprojections[2].norm())){
        PanicButton();
    }
#endif
    } else if (eigenvalues[0]==eigenvalues[1]){
        eigenprojections[2] = eigenvalues[2]/(2.*eigenvalues[2]*eigenvalues[2]*eigenvalues[2]-i1*eigenvalues[2]*eigenvalues[2]+i3) *(aux*aux - (i1-eigenvalues[2])*aux + (i3/eigenvalues[2])*ident);
        eigenprojections[0] = ident-eigenprojections[2];
        eigenprojections[1] = ident-eigenprojections[2];
#ifdef DEBUG_BUILD
    if (std::isnan(eigenprojections[0].norm())||std::isnan(eigenprojections[1].norm())||std::isnan(eigenprojections[2].norm())){
        PanicButton();
    }
#endif
    } else if (eigenvalues[0]==eigenvalues[2]){
        eigenprojections[1] = eigenvalues[1]/(2.*eigenvalues[1]*eigenvalues[1]*eigenvalues[1]-i1*eigenvalues[1]*eigenvalues[1]+i3) *(aux*aux - (i1-eigenvalues[1])*aux + (i3/eigenvalues[1])*ident);
        eigenprojections[0] = ident-eigenprojections[1];
        eigenprojections[2] = ident-eigenprojections[1];
        
#ifdef DEBUG_BUILD
    if (std::isnan(eigenprojections[0].norm())||std::isnan(eigenprojections[1].norm())||std::isnan(eigenprojections[2].norm())){
        PanicButton();
    }
#endif        
    } else if (eigenvalues[1]==eigenvalues[2]) {
        eigenprojections[0] = eigenvalues[0]/(2.*eigenvalues[0]*eigenvalues[0]*eigenvalues[0]-i1*eigenvalues[0]*eigenvalues[0]+i3) *(aux*aux - (i1-eigenvalues[0])*aux + (i3/eigenvalues[0])*ident);
        eigenprojections[2] = ident-eigenprojections[0];
        eigenprojections[1] = ident-eigenprojections[0];
#ifdef DEBUG_BUILD
    if (std::isnan(eigenprojections[0].norm())||std::isnan(eigenprojections[1].norm())||std::isnan(eigenprojections[2].norm())){
        PanicButton();
    }
#endif
    } else {
        eigenprojections[0] = ident;
        eigenprojections[1] = ident;
        eigenprojections[2] = ident;
#ifdef DEBUG_BUILD
    if (std::isnan(eigenprojections[0].norm())||std::isnan(eigenprojections[1].norm())||std::isnan(eigenprojections[2].norm())){
        PanicButton();
    }
#endif
    }

    std::vector<int> newindex(3);
    // Sorting in descending order manually
    if (eigenvalues[0] > eigenvalues[1] && eigenvalues[0] > eigenvalues[2]){
        newindex[0] = 0;
        if (eigenvalues[1] > eigenvalues[2]){
            newindex[1] = 1;
            newindex[2] = 2;
        } else {
            newindex[1] = 2;
            newindex[2] = 1;
        }
    } else if (eigenvalues[1] > eigenvalues[2]){
        newindex[0] = 1;
        if (eigenvalues[0] > eigenvalues[2]){
            newindex[1] = 0;
            newindex[2] = 2;
        } else {
            newindex[1] = 2;
            newindex[2] = 0;
        }
    } else {
        newindex[0] = 2;
        if (eigenvalues[0] > eigenvalues[1]){
            newindex[1] = 0;
            newindex[2] = 1;
        } else {
            newindex[1] = 1;
            newindex[2] = 0;
        }
    }
    auto copyeigval = eigenvalues;
    auto copyeigproj = eigenprojections;
    for (int i = 0; i < 3; i++){
        eigenvalues[i] = copyeigval[newindex[i]];
        eigenprojections[i] = copyeigproj[newindex[i]];
    }

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