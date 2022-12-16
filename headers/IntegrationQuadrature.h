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

template<int DIM, int DEG>
class IntegQuadrature{
public:
    int getNumberOfIntegrationPoints() const {
        return numIntegPoints;
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
    IntegQuadrature(){
        pointWeight.resize(numIntegPoints);
        pointCoord.resize(numIntegPoints,DIM);
        setQuadrature();
    }

private:
    const static int numIntegPoints = (DEG*DEG*(3-DIM)-12*(DIM-2)+DEG*(15*DIM-25))/2;

    ///List of integration points coordinates
    MatrixDouble pointCoord;

    ///List of integration points weights
    VecDouble pointWeight;

};


#endif
