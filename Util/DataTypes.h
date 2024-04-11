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

#include "DenseEigen.h"
#include "SparseEigen.h"
#include "PETScDataTypes.h"

/// @brief Element topologies implemented in the code
enum ElementType {EPoint, EOneD, ETriangle, ETetrahedron, EQuadrilateral, EHexahedron};

#endif 
