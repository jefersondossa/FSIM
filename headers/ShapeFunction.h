//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------QUADRATIC SHAPE FUNCTION----------------------------
//------------------------------------------------------------------------------

#ifndef SHAPEFUNCTION_H
#define SHAPEFUNCTION_H

#include "DataTypes.h"

/// Defines the quadratic shape functions and its derivatives

template<int DIM, int DEG>
class ShapeFunction{
public:
    //Number of element notes    
    static const int numIntpNodes = 4*DIM-2;

    /// Evaluates the shape function value
    /// @param double* Non-dimensional coordinates 
    /// @param double* Shape function values
    void evaluate(VecDouble &xi, VecDouble &phi) const;
    
    /// Evaluates the values of the shape funtion derivatives
    /// @param double* Non-dimensional coordinates 
    /// @param double** Shape function derivatives values
    void evaluateGradient(VecDouble &xi, MatrixDouble &dphi) const;   

    /// Evaluates the values of the shape funtion second derivatives    
    /// @param double** Shape function second derivatives values
    void evaluateHessian(double ***ddphi) const;     

    /// Gets the element nodal coordinates
    /// @param double** Non-dimensional coordinate
    void getCoordinates(double** &coord) const;
};


#endif
