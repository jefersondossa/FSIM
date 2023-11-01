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
#include <vector>
/// Defines the shape functions and its derivatives

class ShapeFunction{
public:
    //Constructor
    ShapeFunction(int dim, int deg){
        fOrder = deg;
    }

    /// Shapes the shape function value
    /// @param double* Non-dimensional coordinates 
    /// @param double* Shape function values
    virtual void Shape(VecDouble &xi, VecDouble &phi) const{
        PanicButton();
    };
    
    /// Shapes the values of the shape funtion derivatives
    /// @param double* Non-dimensional coordinates 
    /// @param double** Shape function derivatives values
    virtual void ShapeGradient(VecDouble &xi, MatrixDouble &dphi) const{
        PanicButton();
    };   

    /// Shapes the values of the shape funtion second derivatives    
    /// @param double** Shape function second derivatives values
    virtual void ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) const{
        PanicButton();
    };     

    /// Gets the element nodal coordinates
    /// @param double** Non-dimensional coordinate
    virtual void getCoordinates(MatrixDouble &coord) const {
        PanicButton();
    };

protected:
    int fOrder;
};


#endif
