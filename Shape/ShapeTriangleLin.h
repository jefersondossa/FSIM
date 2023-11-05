#ifndef SHAPETRIANGLELIN
#define SHAPETRIANGLELIN

#include "DataTypes.h"
#include "IntRuleTriangle.h"

class ShapeTriangleLin{
public:
    // Type definition of integration rule for quadrilateral elements
    typedef IntRuleTriangle LocIntRule;

    const static int Dimension = 2;
    const static int Order = 1;
    const static int NElNodes = 3;
    const static int NCornerNodes = 3;
    
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
