#ifndef __DATATYPES__
#define __DATATYPES__
#include "PanicButton.h"

///\cond
#include <cmath>
#include <stdio.h>
#include <vector>
#include <functional>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
///\endcond

///Eigen data types
using namespace Eigen;

/// @brief Matrix (2d tensor) of doubles
typedef Eigen::MatrixXd MatrixDouble;
/// @brief Matrix (2d tensor) of integers
typedef Eigen::MatrixXi MatrixInt;
/// @brief Sparse matrix of doubles
typedef Eigen::SparseMatrix<double> SparseMat;

/// @brief Vector (1d tensor) of integers
typedef Eigen::Matrix<int64_t, Dynamic, 1> VecInt;
/// @brief Vector (1d tensor) of doubles
typedef Eigen::VectorXd VecDouble;
/// @brief Vector (1d tensor) of booleans
typedef Eigen::Matrix<bool, Dynamic, 1>  VecBool;

/// @brief Element topologies implemented in the code
enum ElementType {EPoint, EOneD, ETriangle, ETetrahedron, EQuadrilateral, EHexahedron};

#endif 
