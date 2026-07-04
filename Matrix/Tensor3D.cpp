#include "Tensor3D.h"
#include "PanicButton.h"


bool compare_head(const REAL& lhs, const REAL& rhs)
{
    return lhs > rhs;
}


Tensor3D::Tensor3D(){
    fData.resize(6);
    fData.setZero();
}

Tensor3D::Tensor3D(const Tensor3D &tensor){
    fData = tensor.fData;
}

Tensor3D::Tensor3D(MatrixDouble &tensor){
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

Tensor3D::Tensor3D(VecDouble &tensor){
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

void Tensor3D::SetData(MatrixDouble &tensor){
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

void Tensor3D::Zero(){
    fData.setZero();
}

Tensor3D Tensor3D::Hydrostatic(){
    Tensor3D fHydrostatic(*this);
    fHydrostatic.Zero();
    REAL val = Trace() / 3.;
    fHydrostatic.fData[XX] = val;
    fHydrostatic.fData[YY] = val;
    fHydrostatic.fData[ZZ] = val;
    return fHydrostatic;
}


Tensor3D Tensor3D::Deviatory(){ 
    Tensor3D fDeviatory = *this;
    REAL val = I1() / 3.;
    fDeviatory.fData[XX] -= val;
    fDeviatory.fData[YY] -= val;
    fDeviatory.fData[ZZ] -= val;
    return fDeviatory;
}

REAL Tensor3D::I1() const{
    return fData[XX] + fData[YY] + fData[ZZ];
}

REAL Tensor3D::I2() const{
    return -(fData[XY] * fData[XY] +
             fData[XZ] * fData[XZ] +
             fData[YZ] * fData[YZ])
           +(fData[XX] * fData[YY] +
             fData[YY] * fData[ZZ] +
             fData[XX] * fData[ZZ]);
}

REAL Tensor3D::I3() const{
    return fData[XX] * fData[YY] * fData[ZZ]
         +(fData[XY] * fData[XZ] * fData[YZ]) * 2.
         -(fData[XX] * fData[YZ] * fData[YZ] +
           fData[YY] * fData[XZ] * fData[XZ] +
           fData[ZZ] * fData[XY] * fData[XY]);
}

REAL Tensor3D::J1(){
    auto dev = Deviatory();
    return dev.fData[XX] + dev.fData[YY] + dev.fData[ZZ];
}

REAL Tensor3D::J2() const{
    REAL i1 = I1();
    REAL j2 = i1*i1 / 3. - I2();
    return j2;
}

REAL Tensor3D::Determinant(){
    return fData[XX] * fData[YY] * fData[ZZ] + fData[XY] * fData[XZ] * fData[YZ]*2. - fData[XZ] * fData[YY] * fData[XZ] -
           fData[XY] * fData[XY] * fData[ZZ] - fData[YZ] * fData[YZ] * fData[XX];
}


REAL Tensor3D::J3(){
    auto fDeviatory = Deviatory();
    REAL det = fDeviatory.fData[XX] * fDeviatory.fData[YY] * fDeviatory.fData[ZZ] +
                 fDeviatory.fData[XY] * fDeviatory.fData[XZ] * fDeviatory.fData[YZ]*2. - 
                 fDeviatory.fData[XZ] * fDeviatory.fData[YY] * fDeviatory.fData[XZ] -
                 fDeviatory.fData[XY] * fDeviatory.fData[XY] * fDeviatory.fData[ZZ] - 
                 fDeviatory.fData[YZ] * fDeviatory.fData[YZ] * fDeviatory.fData[XX];
    return det;
}

REAL Tensor3D::Norm(){
    REAL norm = 0.;
    for (unsigned int i = 0; i < 6; i++) {
        norm += fData[i] * fData[i];
    }
    norm += fData[XY] * fData[XY];
    norm += fData[XZ] * fData[XZ];
    norm += fData[YZ] * fData[YZ];
    return sqrt(norm);
}

REAL Tensor3D::DeviatoryNorm(){
    auto fDeviatory = Deviatory();
    REAL norm = 0.;
    for (unsigned int i = 0; i < 6; i++) {
        norm += fDeviatory.fData[i] * fDeviatory.fData[i];
    }
    norm += fDeviatory.fData[XY] * fDeviatory.fData[XY];
    norm += fDeviatory.fData[XZ] * fDeviatory.fData[XZ];
    norm += fDeviatory.fData[YZ] * fDeviatory.fData[YZ];
    return sqrt(norm);
}

REAL Tensor3D::DoubleContraction(Tensor3D &t){
    return fData[XX]*t.fData[XX] + fData[YY]*t.fData[YY] + fData[ZZ]*t.fData[ZZ] + 
           2.*(fData[XY]*t.fData[XY] + fData[XZ]*t.fData[XZ] + fData[YZ]*t.fData[YZ]);
}

void Tensor3D::Identity(){
    fData.setZero();
    fData[XX] = 1.;
    fData[YY] = 1.;
    fData[ZZ] = 1.;
}

const Tensor3D & Tensor3D::operator*=(const REAL &multipl) {
    int i;
    for (i = 0; i < 6; i++)fData[i] *= multipl;
    return *this;
}

Tensor3D Tensor3D::operator*(const REAL &multipl) const {
    Tensor3D temp(*this);
    return temp *= multipl;
}

const Tensor3D & Tensor3D::operator/=(const REAL &multipl) {
    int i;
    for (i = 0; i < 6; i++)fData[i] /= multipl;
    return *this;
}

Tensor3D Tensor3D::operator/(const REAL &multipl) const {
    Tensor3D temp(*this);
    return temp /= multipl;
}

const Tensor3D & Tensor3D::operator+=(const Tensor3D &sum) {
    int i;
    for (i = 0; i < 6; i++)fData[i] += sum.fData[i];
    return *this;
}

Tensor3D Tensor3D::operator+(const Tensor3D &sum) const {
    Tensor3D temp(*this);
    return temp += sum;
}

const Tensor3D & Tensor3D::operator-=(const Tensor3D &sum) {
    int i;
    for (i = 0; i < 6; i++)fData[i] -= sum.fData[i];
    return *this;
}

Tensor3D Tensor3D::operator-(const Tensor3D &sum) const {
    Tensor3D temp(*this);
    return temp -= sum;
}

Tensor3D Tensor3D::Normalized(){
    Tensor3D temp(*this);
    temp *= 1./Norm();
    return temp;
}

Tensor3D Tensor3D::NormalizedDeviatory(){
    Tensor3D temp = Deviatory();
    temp *= 1./DeviatoryNorm();
    return temp;
}

REAL Tensor3D::Trace() const{
    return fData[XX] + fData[YY] + fData[ZZ];
}

VecDouble Tensor3D::Eigenvalues(){    
    MatrixDouble aux = MatrixForm();
    EigenSolver<MatrixDouble> solver(aux,EigenvaluesOnly);
    VecDouble eigval = solver.eigenvalues().real();

    return eigval;
}

MatrixDouble Tensor3D::Eigenvectors(){    
    MatrixDouble aux = MatrixForm();
    
    EigenSolver<MatrixDouble> solver(aux);
    MatrixDouble eigvec = solver.eigenvectors().real();

    return eigvec;
}

void Tensor3D::SpectralDecomposition(VecDouble &eigenvalues, std::vector<MatrixDouble> &eigenprojections){ 
    eigenvalues.resize(3); 
    eigenprojections.resize(3);

    //A.4 Souza Neto
    REAL i1 = I1();
    REAL i2 = I2();
    REAL i3 = I3();
    REAL R = (-2.*i1*i1*i1 + 9.*i1*i2 - 27.*i3) / 54.;
    REAL Q = (i1*i1 - 3.*i2)/9.;
    REAL sqQ = sqrt(Q);
    REAL valcos = R/(Q*sqQ);
    if (valcos < -1.) valcos = -1.;
    if (valcos >  1.) valcos = 1.;
    REAL theta = acos(valcos);
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

Tensor3D Tensor3D::Multiply(MatrixDouble &mat){
    Tensor3D temp(*this);
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

MatrixDouble Tensor3D::MatrixForm(){
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

MatrixDouble Tensor3D::TensorProduct(Tensor3D &tensor){
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