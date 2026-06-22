#include "ArlequinMatRedEigen.h"


ArlequinMatRedEigen::ArlequinMatRedEigen(int dim0, int dim1, int dimLagrange) : MatrixType(dim0+dim1+dimLagrange,dim0+dim1+dimLagrange){
    fDim0 = dim0;
    fDim1 = dim1;
    fDimLagrange = dimLagrange;
    fK0.resize(dim0,dim0);
    fK1.resize(dim1,dim1);
    fL0.resize(dim0,dimLagrange);
    fL1.resize(dim1,dimLagrange);
    fL0T.resize(dimLagrange,dim0);
    fL1T.resize(dimLagrange,dim1);
    fE.resize(dimLagrange,dimLagrange);
    fU0.resize(dim0, 1);
    fU1.resize(dim1, 1);
    fLambda.resize(dimLagrange, 1);
    fF0.resize(dim0, 1);
    fF1.resize(dim1, 1);
    fG.resize(dimLagrange, 1);
    ZeroMatrix();
    ZeroRhs();
    ZeroSolution();
}

ArlequinMatRedEigen::~ArlequinMatRedEigen()
{
    fK0.resize(0,0);
    fK1.resize(0,0);
    fL0.resize(0,0);
    fL1.resize(0,0);
    fL0T.resize(0,0);
    fL1T.resize(0,0);
    fE.resize(0,0);
    fU0.resize(0,0);
    fU1.resize(0,0);
    fLambda.resize(0,0);
    fF0.resize(0,0);
    fF1.resize(0,0);
    fG.resize(0,0);
}

void ArlequinMatRedEigen::ZeroMatrix(){
    fK0.setZero();
    fK1.setZero();
    fL0.setZero();
    fL1.setZero();
    fL0T.setZero();
    fL1T.setZero();
    fE.setZero();
}

void ArlequinMatRedEigen::ZeroSolution(){
    fU0.setZero();
    fU1.setZero();
    fLambda.setZero();
}

void ArlequinMatRedEigen::ZeroRhs(){
    fF0.setZero();
    fF1.setZero();
    fG.setZero();
}

void ArlequinMatRedEigen::AddValueMatrix(int &row, int &col, double &val) {
    if (row < fDim0 && col < fDim0){// Belongs to K0
        fK0.coeffRef(row, col) += val;
        return;
    } else if (row < fDim0 && col > fDim0+fDim1-1){ // Belongs to L0
        fL0.coeffRef(row, col-fDim0-fDim1) += val;
        return;
    } else if (col < fDim0 && row > fDim0+fDim1-1){// Belongs to L0T
        fL0T.coeffRef(row-fDim0-fDim1,col) += val;
        return;
    } else if (row > fDim0-1 && row < fDim0+fDim1 && col > fDim0-1 && col < fDim0+fDim1){ // Belongs to K1
        fK1.coeffRef(row-fDim0, col-fDim0) += val;
        return;
    } else if (row > fDim0-1 && row < fDim0+fDim1 && col > fDim0+fDim1-1){ // Belongs to L1
        fL1.coeffRef(row-fDim0, col-fDim0-fDim1) += val;
        return;
    } else if (col > fDim0-1 && col < fDim0+fDim1 && row > fDim0+fDim1-1){ // Belongs to L1T
        fL1T.coeffRef(row-fDim0-fDim1,col-fDim0) += val;
        return;
    } else if (row > fDim0+fDim1-1 && col > fDim0+fDim1-1){ // Belongs to E
        fE.coeffRef(row-fDim0-fDim1, col-fDim0-fDim1) += val;
        return;
    } else {
        PanicButton();
    }
};
void ArlequinMatRedEigen::PutValueMatrix(int &row, int &col, double &val) {
    if (row < fDim0 && col < fDim0){// Belongs to K0
        fK0.coeffRef(row, col) = val;
        return;
    } else if (row < fDim0 && col > fDim0+fDim1){ // Belongs to L0
        fL0.coeffRef(row, col-fDim0-fDim1) = val;
        return;
    } else if (row >= fDim0 && row < fDim0+fDim1 && col >= fDim0 && col < fDim0+fDim1){ // Belongs to K1
        fK1.coeffRef(row-fDim0, col-fDim0) = val;
        return;
    } else if (row >= fDim0 && row < fDim0+fDim1 && col > fDim0+fDim1){ // Belongs to L1
        fL1.coeffRef(row-fDim0, col-fDim0-fDim1) = val;
        return;
    } else if (row > fDim0+fDim1 && col > fDim0+fDim1){ // Belongs to E
        fE.coeffRef(row-fDim0-fDim1, col-fDim0-fDim1) = val;
        return;
    } else {
        PanicButton();
    }
};
double ArlequinMatRedEigen::GetValueMatrix(int &row, int &col) {
    if (row < fDim0 && col < fDim0){// Belongs to K0
        return fK0.coeffRef(row, col);
    } else if (row < fDim0 && col > fDim0+fDim1){ // Belongs to L0
        return fL0.coeffRef(row, col-fDim0-fDim1);
    } else if (row >= fDim0 && row < fDim0+fDim1 && col >= fDim0 && col < fDim0+fDim1){ // Belongs to K1
        return fK1.coeffRef(row-fDim0, col-fDim0);
    } else if (row >= fDim0 && row < fDim0+fDim1 && col > fDim0+fDim1){ // Belongs to L1
        return fL1.coeffRef(row-fDim0, col-fDim0-fDim1);
    } else if (row > fDim0+fDim1 && col > fDim0+fDim1){ // Belongs to E
        return fE.coeffRef(row-fDim0-fDim1, col-fDim0-fDim1);
    } else {
        PanicButton();
    }
};

void ArlequinMatRedEigen::AddValueRhs(int &row, double &val) {
    if (row < fDim0){// Belongs to F0
        fF0(row,0) += val;
        return;
    } else if (row >= fDim0 && row < fDim0+fDim1){ // Belongs to F1
        fF1(row-fDim0, 0) += val;
        return;
    } else if (row >= fDim0+fDim1){ // Belongs to G
        fG(row-fDim0-fDim1, 0) += val;
        return;
    } else {
        PanicButton();
    }   
};
void ArlequinMatRedEigen::PutValueRhs(int &row, double &val) {
    if (row < fDim0){// Belongs to F0
        fF0(row,0) = val;
        return;
    } else if (row >= fDim0 && row < fDim0+fDim1){ // Belongs to F1
        fF1(row-fDim0, 0) = val;
        return;
    } else if (row >= fDim0+fDim1){ // Belongs to G
        fG(row-fDim0-fDim1, 0) = val;
        return;
    } else {
        PanicButton();
    }   
};
double ArlequinMatRedEigen::GetValueRhs(int &row) {
    if (row < fDim0){// Belongs to F0
        return fF0(row,0);
    } else if (row >= fDim0 && row < fDim0+fDim1){ // Belongs to F1
        return fF1(row-fDim0, 0);
    } else if (row >= fDim0+fDim1){ // Belongs to G
        return fG(row-fDim0-fDim1, 0);
    } else {
        PanicButton();
    }   
};

void ArlequinMatRedEigen::AddValueSolution(int &row, double &val) {
    if (row < fDim0){// Belongs to U0
        fU0(row,0) += val;
        return;
    } else if (row >= fDim0 && row < fDim0+fDim1){ // Belongs to U1
        fU1(row-fDim0, 0) += val;
        return;
    } else if (row >= fDim0+fDim1){ // Belongs to Lambda
        fLambda(row-fDim0-fDim1, 0) += val;
        return;
    } else {
        PanicButton();
    }   
};
void ArlequinMatRedEigen::PutValueSolution(int &row, double &val) {
    if (row < fDim0){// Belongs to U0
        fU0(row,0) = val;
        return;
    } else if (row >= fDim0 && row < fDim0+fDim1){ // Belongs to U1
        fU1(row-fDim0, 0) = val;
        return;
    } else if (row >= fDim0+fDim1){ // Belongs to Lambda
        fLambda(row-fDim0-fDim1, 0) = val;
        return;
    } else {
        PanicButton();
    }   
};
double ArlequinMatRedEigen::GetValueSolution(int &row) {
    if (row < fDim0){// Belongs to U0
        return fU0(row,0);
    } else if (row >= fDim0 && row < fDim0+fDim1){ // Belongs to U1
        return fU1(row-fDim0, 0);
    } else if (row >= fDim0+fDim1){ // Belongs to Lambda
        return fLambda(row-fDim0-fDim1, 0);
    } else {
        PanicButton();
    }   
};

double ArlequinMatRedEigen::SolutionNorm(){
    return fU0.norm()+fU1.norm();
}

void ArlequinMatRedEigen::PrintMatrix(){
    // std::cout << "Global Matrix = \n"<< fMatrix << std::endl;
    std::cout << "K0 = {\n";
    for (int i = 0; i < fK0.rows(); i++){
        std::cout << "{";
        for (int j = 0; j < fK0.cols(); j++){
            std::cout << fK0.coeffRef(i,j);
            if (j<fK0.cols()-1) std::cout << ",";
        }
        std::cout << "}";
        if (i<fK0.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 

    std::cout << "K1 = {\n";
    for (int i = 0; i < fK1.rows(); i++){
        std::cout << "{";
        for (int j = 0; j < fK1.cols(); j++){
            std::cout << fK1.coeffRef(i,j);
            if (j<fK1.cols()-1) std::cout << ",";
        }
        std::cout << "}";
        if (i<fK1.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n";
    
    std::cout << "L0 = {\n";
    for (int i = 0; i < fL0.rows(); i++){
        std::cout << "{";
        for (int j = 0; j < fL0.cols(); j++){
            std::cout << fL0.coeffRef(i,j);
            if (j<fL0.cols()-1) std::cout << ",";
        }
        std::cout << "}";
        if (i<fL0.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 

    std::cout << "L1 = {\n";
    for (int i = 0; i < fL1.rows(); i++){
        std::cout << "{";
        for (int j = 0; j < fL1.cols(); j++){
            std::cout << fL1.coeffRef(i,j);
            if (j<fL1.cols()-1) std::cout << ",";
        }
        std::cout << "}";
        if (i<fL1.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 

    std::cout << "E = {\n";
    for (int i = 0; i < fE.rows(); i++){
        std::cout << "{";
        for (int j = 0; j < fE.cols(); j++){
            std::cout << fE.coeffRef(i,j);
            if (j<fE.cols()-1) std::cout << ",";
        }
        std::cout << "}";
        if (i<fE.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 
}

void ArlequinMatRedEigen::PrintRhs(){
    // std::cout << "Global RHS = \n"<< fRhs << std::endl;
    std::cout << "F0 = {\n";
    for (int i = 0; i < fF0.rows(); i++){
        std::cout << "{" << fF0(i,0) << "}";
        if (i<fF0.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n";

    std::cout << "F1 = {\n";
    for (int i = 0; i < fF1.rows(); i++){
        std::cout << "{" << fF1(i,0) << "}";
        if (i<fF1.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 

    std::cout << "G = {\n";
    for (int i = 0; i < fG.rows(); i++){
        std::cout << "{" << fG(i,0) << "}";
        if (i<fG.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 
    
}

void ArlequinMatRedEigen::PrintSolution(){
    // std::cout << "Solution = \n"<< fSolution << std::endl;
    std::cout << "U0 = {\n";
    for (int i = 0; i < fU0.rows(); i++){
        std::cout << "{" << fU0(i,0) << "}";
        if (i<fU0.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n";

    std::cout << "U1 = {\n";
    for (int i = 0; i < fU1.rows(); i++){
        std::cout << "{" << fU1(i,0) << "}";
        if (i<fU1.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n";

    std::cout << "Lambda = {\n";
    for (int i = 0; i < fLambda.rows(); i++){
        std::cout << "{" << fLambda(i,0) << "}";
        if (i<fLambda.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 
}