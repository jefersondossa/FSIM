#ifndef DENSE_EIGEN
#define DENSE_EIGEN

#include <cmath>
#include <stdio.h>
#include <vector>
#include <functional>
#include <Eigen/Dense>
#include <iostream>
#include "config.h"

using namespace Eigen;

#ifdef REALfloat
typedef float REAL;
#endif // REALfloat
#ifdef REALdouble
typedef double REAL; //This is the default configuration
#endif // REALdouble
#ifdef REALlongdouble
typedef long double REAL;
#endif // REALlongdouble

/// @brief Matrix (2d tensor) of doubles
typedef Eigen::Matrix<REAL, -1, -1> MatrixDouble;
typedef Eigen::Matrix<REAL, 2, 2> Matrix2R;
/// @brief Matrix (2d tensor) of integers
typedef Eigen::MatrixXi MatrixInt;

/// @brief Vector (1d tensor) of integers
typedef Eigen::Matrix<int, Dynamic, 1> VecInt;
/// @brief Vector (1d tensor) of doubles
typedef Eigen::Matrix<REAL, -1, 1> VecDouble;
/// @brief Vector (1d tensor) of 3D vectors
typedef Eigen::Matrix<REAL, 3, 1> Vector3R;
/// @brief Vector (1d tensor) of booleans
typedef Eigen::Matrix<bool, Dynamic, 1>  VecBool;

void PrintMathematica(MatrixDouble &mat, std::string name);
void PrintMathematica(VecDouble &vec, std::string name);

#endif