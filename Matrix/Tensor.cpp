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
    MatrixDouble aux(3,3);
    aux(0,0) = fData[XX];
    aux(1,1) = fData[YY];
    aux(2,2) = fData[ZZ];
    aux(0,1) = aux(1,0) = fData[XY];
    aux(0,2) = aux(2,0) = fData[XZ];
    aux(1,2) = aux(2,1) = fData[YZ];

    EigenSolver<MatrixDouble> solver(aux,EigenvaluesOnly);
    VecDouble eigval = solver.eigenvalues().real();

    std::vector<double> auxvec(3);
    auxvec[0] = eigval[0];
    auxvec[1] = eigval[1];
    auxvec[2] = eigval[2];
    //Sort in ascending order
    std::sort(auxvec.begin(), auxvec.end(), &compare_head);

    for (int64_t i = 0; i < eigval.rows(); ++i)
        eigval[i] = auxvec[i];


    //A.4 Souza Neto
    double i1 = I1();
    double i2 = I2();
    double i3 = I3();
    double R = (-2.*i1*i1*i1 + 9.*i1*i2 - 27.*i3) / 54.;
    double Q = (i1*i1 - 3.*i2)/9.;
    double sqQ = sqrt(Q);
    double theta = acos(R/(Q*sqQ));
    double val1 = -2. * sqQ * cos(theta/3.) + i1/3.; 
    double val2 = -2. * sqQ * cos((theta+2.*M_PI)/3.) + i1/3.; 
    double val3 = -2. * sqQ * cos((theta-2.*M_PI)/3.) + i1/3.; 

    return eigval;
}

void Tensor::SpectralDecomposition(VecDouble &eigenvalues, MatrixDouble &eigenvectors){ 
    eigenvalues.resize(3); 
    eigenvectors.resize(3,3);

    MatrixDouble aux(3,3);
    aux = MatrixForm();
    // aux(0,0) = fData[XX];
    // aux(1,1) = fData[YY];
    // aux(2,2) = fData[ZZ];
    // aux(0,1) = aux(1,0) = fData[XY];
    // aux(0,2) = aux(2,0) = fData[XZ];
    // aux(1,2) = aux(2,1) = fData[YZ];

    EigenSolver<MatrixDouble> solver(aux);
    eigenvalues = solver.eigenvalues().real();
    eigenvectors = solver.eigenvectors().real();

    std::vector<int> newindex(3);
    //Sorting in descending order manually
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
    auto copyeigvec = eigenvectors;
    for (int i = 0; i < 3; i++){
        eigenvalues[i] = copyeigval[newindex[i]];
        for (int j = 0; j < 3; j++){
            eigenvectors(i,j) = copyeigvec(newindex[i],newindex[j]);   
        }
    }
    //A.4 Souza Neto
    double i1 = I1();
    double i2 = I2();
    double i3 = I3();
    double R = (-2.*i1*i1*i1 + 9.*i1*i2 - 27.*i3) / 54.;
    double Q = (i1*i1 - 3.*i2)/9.;
    double sqQ = sqrt(Q);
    double theta = acos(R/(Q*sqQ));
    VecDouble xi(3);
    xi[0] = -2. * sqQ * cos(theta/3.) + i1/3.; 
    xi[1] = -2. * sqQ * cos((theta+2.*M_PI)/3.) + i1/3.; 
    xi[2] = -2. * sqQ * cos((theta-2.*M_PI)/3.) + i1/3.;

    // If all eigenvalues are different
    MatrixDouble ident(3,3);
    ident.setIdentity();
    std::vector<MatrixDouble> eigenprojections(3);
    if (fabs(xi[0]-xi[1]) > 1.e-10 && fabs(xi[0]-xi[2]) > 1.e-10 && fabs(xi[1]-xi[2]) > 1.e-10){
        for (size_t i = 0; i < 3; i++){
            eigenprojections[i] = xi[i]/(2.*xi[i]*xi[i]*xi[i]-i1*xi[i]*xi[i]+i3) *(aux*aux - (i1-xi[i])*aux + (i3/xi[i])*ident);
        } // 
    } else if (fabs(xi[0]-xi[1]) > 1.e-10){
        eigenprojections[0] = xi[0]/(2.*xi[0]*xi[0]*xi[0]-i1*xi[0]*xi[0]+i3) *(aux*aux - (i1-xi[0])*aux + (i3/xi[0])*ident);
        eigenprojections[1] = eigenprojections[0] - ident;
        eigenprojections[2] = xi[2]/(2.*xi[2]*xi[2]*xi[2]-i1*xi[2]*xi[2]+i3) *(aux*aux - (i1-xi[2])*aux + (i3/xi[2])*ident);
    } else if (fabs(xi[0]-xi[2]) > 1.e-10){
        eigenprojections[0] = xi[0]/(2.*xi[0]*xi[0]*xi[0]-i1*xi[0]*xi[0]+i3) *(aux*aux - (i1-xi[0])*aux + (i3/xi[0])*ident);
        eigenprojections[2] = eigenprojections[0] - ident;
        eigenprojections[1] = xi[1]/(2.*xi[1]*xi[1]*xi[1]-i1*xi[1]*xi[1]+i3) *(aux*aux - (i1-xi[1])*aux + (i3/xi[1])*ident);
    } else if (fabs(xi[1]-xi[2]) > 1.e-10) {
        eigenprojections[0] = xi[0]/(2.*xi[0]*xi[0]*xi[0]-i1*xi[0]*xi[0]+i3) *(aux*aux - (i1-xi[0])*aux + (i3/xi[0])*ident);
        eigenprojections[2] = xi[2]/(2.*xi[2]*xi[2]*xi[2]-i1*xi[2]*xi[2]+i3) *(aux*aux - (i1-xi[2])*aux + (i3/xi[2])*ident);
        eigenprojections[1] = eigenprojections[2] - ident;
    } else {
        eigenprojections[0] = ident;
        eigenprojections[1] = ident;
        eigenprojections[2] = ident;
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