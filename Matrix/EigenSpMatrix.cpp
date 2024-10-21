#include "EigenSpMatrix.h"


EigenSpMatrix::EigenSpMatrix(int64_t rows, int64_t cols) : MatrixType(rows,cols){
    
    fMatrix.resize(rows, cols);
    fMatrix.reserve(VectorXi::Constant(cols, 128));
    fRhs.resize(rows, 1);
    fSolution.resize(rows, 1);
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

void EigenSpMatrix::AddValueMatrix(int64_t &row, int64_t &col, double &val) {
    fMatrix.coeffRef(row, col) += val;
};
void EigenSpMatrix::PutValueMatrix(int64_t &row, int64_t &col, double &val) {
    fMatrix.coeffRef(row, col) = val;
};
double EigenSpMatrix::GetValueMatrix(int64_t &row, int64_t &col) {
    return fMatrix.coeffRef(row, col);
};

void EigenSpMatrix::AddValueRhs(int64_t &row, double &val) {
    fRhs(row,0) += val;
};
void EigenSpMatrix::PutValueRhs(int64_t &row, double &val) {
    fRhs(row,0) = val;
};
double EigenSpMatrix::GetValueRhs(int64_t &row) {
    return fRhs(row,0);
};

void EigenSpMatrix::AddValueSolution(int64_t &row, double &val) {
    fSolution(row,0) += val;
};
void EigenSpMatrix::PutValueSolution(int64_t &row, double &val) {
    fSolution(row,0) = val;
};
double EigenSpMatrix::GetValueSolution(int64_t &row) {
    return fSolution(row,0);
};

double EigenSpMatrix::SolutionNorm(){
    return fSolution.norm();
}

void EigenSpMatrix::PrintMatrix(){
    // std::cout << "Global Matrix = \n"<< fMatrix << std::endl;
    std::cout << "GlobalMatrix = {\n";
    for (int64_t i = 0; i < fMatrix.rows(); i++){
        std::cout << "{";
        for (int64_t j = 0; j < fMatrix.cols(); j++){
            std::cout << fMatrix.coeffRef(i,j);
            if (j<fMatrix.cols()-1) std::cout << ",";
        }
        std::cout << "}";
        if (i<fMatrix.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 
}

void EigenSpMatrix::PrintRhs(){
    // std::cout << "Global RHS = \n"<< fRhs << std::endl;
    std::cout << "RHS = {\n";
    for (int64_t i = 0; i < fRhs.rows(); i++){
        std::cout << "{" << fRhs(i,0) << "}";
        if (i<fRhs.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 
    
}

void EigenSpMatrix::PrintSolution(){
    // std::cout << "Solution = \n"<< fSolution << std::endl;
    std::cout << "Solution = {\n";
    for (int64_t i = 0; i < fSolution.rows(); i++){
        std::cout << "{" << fSolution(i,0) << "}";
        if (i<fSolution.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 
}