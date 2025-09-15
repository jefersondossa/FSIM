#ifndef DENSE_EIGEN
#define DENSE_EIGEN

#include <cmath>
#include <stdio.h>
#include <vector>
#include <functional>
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

/// @brief Matrix (2d tensor) of doubles
typedef Eigen::MatrixXd MatrixDouble;
/// @brief Matrix (2d tensor) of integers
typedef Eigen::MatrixXi MatrixInt;

/// @brief Vector (1d tensor) of integers
typedef Eigen::Matrix<int64_t, Dynamic, 1> VecInt;
/// @brief Vector (1d tensor) of doubles
typedef Eigen::VectorXd VecDouble;
/// @brief Vector (1d tensor) of booleans
typedef Eigen::Matrix<bool, Dynamic, 1>  VecBool;

void PrintMathematica(MatrixDouble &mat, std::string name);
void PrintMathematica(VecDouble &vec, std::string name);

#endif