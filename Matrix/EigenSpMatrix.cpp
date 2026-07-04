#include "EigenSpMatrix.h"


EigenSpMatrix::EigenSpMatrix(int rows, int cols) : MatrixType(rows,cols){
    
    fMatrix.resize(rows, cols);
    fMatrix.reserve(VectorXi::Constant(cols, 128));
    fRhs.resize(rows, 1);
    fSolution.resize(rows, 1);
    fRhs.setZero();
    fSolution.setZero();
}

EigenSpMatrix::~EigenSpMatrix()
{
    fMatrix.resize(0,0);
    fRhs.resize(0,0);
    fSolution.resize(0,0);
}

void EigenSpMatrix::ZeroMatrix(){
    fMatrix.setZero();
}

void EigenSpMatrix::ZeroSolution(){
    fSolution.setZero();
}

void EigenSpMatrix::ZeroRhs(){
    fRhs.setZero();
}

void EigenSpMatrix::AddValueMatrix(int &row, int &col, REAL &val) {
    fMatrix.coeffRef(row, col) += val;
};
void EigenSpMatrix::PutValueMatrix(int &row, int &col, REAL &val) {
    fMatrix.coeffRef(row, col) = val;
};
REAL EigenSpMatrix::GetValueMatrix(int &row, int &col) {
    return fMatrix.coeffRef(row, col);
};

void EigenSpMatrix::AddValueRhs(int &row, REAL &val) {
    fRhs(row,0) += val;
};
void EigenSpMatrix::PutValueRhs(int &row, REAL &val) {
    fRhs(row,0) = val;
};
REAL EigenSpMatrix::GetValueRhs(int &row) {
    return fRhs(row,0);
};

void EigenSpMatrix::AddValueSolution(int &row, REAL &val) {
    fSolution(row,0) += val;
};
void EigenSpMatrix::PutValueSolution(int &row, REAL &val) {
    fSolution(row,0) = val;
};
REAL EigenSpMatrix::GetValueSolution(int &row) {
    return fSolution(row,0);
};

REAL EigenSpMatrix::SolutionNorm(){
    return fSolution.norm();
}

void EigenSpMatrix::PrintMatrix(){
    // std::cout << "Global Matrix = \n"<< fMatrix << std::endl;
    std::cout << "GlobalMatrix = {\n";
    for (int i = 0; i < fMatrix.rows(); i++){
        std::cout << "{";
        for (int j = 0; j < fMatrix.cols(); j++){
            if (fabs(fMatrix.coeffRef(i,j))<1.e-13) {
                std::cout << 0.;
            } else {
                std::cout << fMatrix.coeffRef(i,j);
            }
            if (j<fMatrix.cols()-1) std::cout << ",";
        }
        std::cout << "}";
        if (i<fMatrix.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 

    // std::cout << "Matrix Diagonal = {\n";
    // for (int i = 0; i < fMatrix.rows(); i++){
    //     std::cout << "{";
    //     if (fabs(fMatrix.coeffRef(i,i))<1.e-13) {
    //         std::cout << 0.;
    //     } else {
    //         std::cout << fMatrix.coeffRef(i,i);
    //     }
    //     std::cout << "}";
    //     if (i<fMatrix.rows()-1) std::cout << ",\n";
    // }
    // std::cout << "};\n"; 
}

void EigenSpMatrix::PrintRhs(){
    // std::cout << "Global RHS = \n"<< fRhs << std::endl;
    std::cout << "RHS = {\n";
    for (int i = 0; i < fRhs.rows(); i++){
        std::cout << "{";
        if (fabs(fRhs(i,0))<1.e-13) {
            std::cout << 0.;
        } else {
            std::cout << fRhs(i,0);
        }
        std::cout << "}";
        if (i<fRhs.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 
    
}

void EigenSpMatrix::PrintSolution(){
    // std::cout << "Solution = \n"<< fSolution << std::endl;
    std::cout << "Solution = {\n";
    for (int i = 0; i < fSolution.rows(); i++){
        std::cout << "{";
        if (fabs(fSolution(i,0))<1.e-13) {
            std::cout << 0.;
        } else {
            std::cout << fSolution(i,0);
        }
        std::cout << "}";
        if (i<fSolution.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 
}