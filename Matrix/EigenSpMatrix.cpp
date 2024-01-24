#include "EigenSpMatrix.h"


EigenSpMatrix::EigenSpMatrix(int64_t rows, int64_t cols) : MatrixType(rows,cols){
    
    fMatrix.resize(rows, cols);
    fRhs.resize(rows, 1);
    fSolution.resize(rows, 1);
    ZeroMatrix();
    ZeroRhs();
    ZeroSolution();
}

EigenSpMatrix::~EigenSpMatrix()
{

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
    std::cout << "Global Matrix = \n"<< fMatrix << std::endl;
}

void EigenSpMatrix::PrintRhs(){
    std::cout << "Global RHS = \n"<< fRhs << std::endl;
}

void EigenSpMatrix::PrintSolution(){
    std::cout << "Solution = \n"<< fSolution << std::endl;
}