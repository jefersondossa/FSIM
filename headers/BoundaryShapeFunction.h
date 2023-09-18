//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2020 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------BOUNDARY SHAPE FUNCTIONS----------------------------
//------------------------------------------------------------------------------

#ifndef BOUND_SHAPEFUNCTION_H
#define BOUND_SHAPEFUNCTION_H

#include <iterator>
#include "DataTypes.h"

/// Defines the fluid boundary shape functions
namespace BoundShapeFunction {
    /// Conputes the boundary shape function value
    /// @param double Adimensional coordinates 
    /// @param double* boundary shape function
    /// @param double* boundary shape function derivatives
    void getShapeFunction(int DIM, int DEG, VecDouble &Xsi, VecDouble &phi_, MatrixDouble &dphi_);
};


#endif
