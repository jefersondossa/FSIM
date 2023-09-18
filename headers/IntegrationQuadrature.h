//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------QUADRATURE POINTS-------------------------------
//------------------------------------------------------------------------------

#ifndef INTEG_QUADRATURE_H
#define INTEG_QUADRATURE_H

#include "ShapeFunction.h"
#include <iterator>
#include <math.h>
#include "DataTypes.h"

/// Defines the domain integration Hammer quadrature

class IntegQuadrature{
public:
    int getNumberOfIntegrationPoints() const {
        return pointWeight.size();
    }

    /// Returns the integration point coordinate
    /// @param int integration point index @param int adimensional direction
    /// @return integration point adimensional coordinates
    double PointList(int i, int j) const {return pointCoord(i,j);};

    /// Retuns the integration point weight
    /// @param int integration point index @return integration point weight
    double WeightList(int i) const {return pointWeight[i];};

    /// Sets the domain quadrature and store the values of the integration points and weights 
    void setQuadrature();

    /// Constructor of the domain integration quadrature
    IntegQuadrature(int dim, int deg){
        int numIntegPoints = (deg*deg*(3-dim)-12*(dim-2)+deg*(15*dim-25))/2;
        fOrder = deg;
        pointWeight.resize(numIntegPoints);
        pointCoord.resize(numIntegPoints,dim);
        setQuadrature();
    }

private:
    int fOrder = 0;

    ///List of integration points coordinates
    MatrixDouble pointCoord;

    ///List of integration points weights
    VecDouble pointWeight;

};


#endif
