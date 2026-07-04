#ifndef SPARSE_EIGEN_H
#define SPARSE_EIGEN_H

#include <Eigen/SparseCore>
#include "DenseEigen.h"

using namespace Eigen;

/// @brief Sparse matrix of doubles
typedef Eigen::SparseMatrix<REAL> SparseMat;

#endif