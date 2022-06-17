//
//  IntRule.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __DATATYPES__
#define __DATATYPES__

///\cond
#include <cmath>
#include <stdio.h>
#include <vector>
#include <functional>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
///\endcond

using namespace Eigen;

typedef Eigen::MatrixXd MatrixDouble;
typedef Eigen::MatrixXi MatrixInt;
typedef Eigen::SparseMatrix<double> SparseMat;

//typedef Eigen::VectorXi VecInt;
typedef Eigen::Matrix<int64_t, Dynamic, 1> VecInt;
typedef Eigen::VectorXd VecDouble;
typedef Eigen::Matrix<bool, Dynamic, 1>  VecBool;

#endif /* defined(__FemSC__DATATYPES__) */
