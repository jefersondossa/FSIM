#ifndef SHAPETRIANGLEQUA
#define SHAPETRIANGLEQUA

#include "DataTypes.h"
#include "IntRuleTriangle.h"

class ShapeTriangleQua{
public:
    // Type definition of integration rule for quadrilateral elements
    typedef IntRuleTriangle LocIntRule;

    const static int Dimension = 2;
    const static int Order = 2;
    const static int NElNodes = 6;

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
