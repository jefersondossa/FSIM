#ifndef ShapeOneDCub_H
#define ShapeOneDCub_H

#include "DataTypes.h"
#include "IntRule1d.h"

class ShapeOneDCub{
public:
    // Type definition for integration rule 1D
    typedef IntRule1d LocIntRule;

    const static int Dimension = 1;
    const static int Order = 3;
    const static int NElNodes = 4;
    const static int NCornerNodes = 2;
    
    /// Shapes the shape function value
    /// @param double* Non-dimensional coordinates 
    /// @param double* Shape function values
    static void Shape(VecDouble &xi, VecDouble &phi);
    
    /// Shapes the values of the shape funtion derivatives
    /// @param double* Non-dimensional coordinates 
    /// @param double** Shape function derivatives values
    static void ShapeGradient(VecDouble &xi, MatrixDouble &dphi);

    /// Shapes the values of the shape funtion second derivatives    
    /// @param double** Shape function second derivatives values
    static void ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi);

    /// Gets the element nodal coordinates
    /// @param double** Non-dimensional coordinate
    static void getCoordinates(MatrixDouble &coord);
   
};

#endif
