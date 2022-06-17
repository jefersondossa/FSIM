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

#ifndef INTEG_QUADRATURE11_H
#define INTEG_QUADRATURE11_H

#include "QuadraticShapeFunction.h"
#include <math.h>
#include <iterator>

/// Defines a special quadrature rule (equal to the normal quadradure rule) used for the integration of problems with coincident meshes

template<int DIM, int DEG>
class IntegQuadratureSpecial{
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
    IntegQuadratureSpecial(){
        pointWeight.resize(numIntegPoints);
        pointCoord.resize(numIntegPoints,DIM);
        setQuadrature();
    }
  
    //Interpolate variables
    double interpolateQuadraticVariable(VecDouble &nValues, int point) const ;

private:
    const static int numIntegPoints = -5*DIM-8*DEG+6*DIM*DEG+9;

    ///List of integration points coordinates
    MatrixDouble pointCoord;

    ///List of integration points weights
    VecDouble pointWeight;

    //Defines shape functions
    QuadShapeFunction<DIM,DEG> shapeQuad;

};


#endif
