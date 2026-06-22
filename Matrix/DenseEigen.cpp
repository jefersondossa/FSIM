#include "DenseEigen.h"

void PrintMathematica(MatrixDouble &mat, std::string name = "A"){
    std::cout << name << " = {";
    for (int i = 0; i < mat.rows(); i++){
        std::cout << "{";
        for (int j = 0; j < mat.cols(); j++){
            if (fabs(mat(i,j))<1.e-13) {
                std::cout << 0.;
            } else {
                std::cout << mat(i,j);
            }
            if (j<mat.cols()-1) std::cout << ",";
        }
        std::cout << "}";
        if (i<mat.rows()-1) std::cout << ",\n";
    }
    std::cout << "};\n"; 
};

void PrintMathematica(VecDouble &vec, std::string name = "b"){
    std::cout << name << " = {";
    for (int i = 0; i < vec.size(); i++) {
        if (fabs(vec(i))<1.e-13) {
            std::cout << 0.;
        } else {
            std::cout << vec(i);
        }
        if (i<vec.size()-1) std::cout << ",";
    }
    std::cout << "};\n";
}