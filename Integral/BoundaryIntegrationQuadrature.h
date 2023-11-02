//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2020 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------------BOUNDARY QUADRATURE POINTS---------------------------
//------------------------------------------------------------------------------

#ifndef BOUND_INTEG_QUADRATURE_H
#define BOUND_INTEG_QUADRATURE_H

#include "BoundaryShapeFunction.h"
#include <math.h>
#include "DataTypes.h"

/// Defines the quadrature rule for the boundary integration
class BoundaryIntegQuadrature{
private:
    double pi = M_PI;  //Pi

public:
    int NPoints() {
        return pointWeight.size();
    }

    /// Returns the integration point coordinate
    /// @param int integration point index @param int adimensional direction
    /// @return integration point adimensional coordinates
    double PointList(int i, int j) {return pointCoord(i,j);}; 

    /// Retuns the integration point weight
    /// @param int integration point index @return integration point weight
    double WeightList(int i) {return pointWeight[i];};

    /// Sets the integration quadrature
    void setQuadrature();

    /// Boundary integration quadrature constructor
    BoundaryIntegQuadrature(int DIM, int DEG){
        fOrder = DEG;
        pointWeight.resize(-3*DIM-8*DEG+4*DIM*DEG+8);
        pointCoord.resize(-3*DIM-8*DEG+4*DIM*DEG+8,DIM-1);
        setQuadrature();
    }


private:
    ///List of integration points coordinates
    MatrixDouble pointCoord;

    ///List of integration points weights
    VecDouble pointWeight;

    int fOrder;
};


#endif
