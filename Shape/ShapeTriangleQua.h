#ifndef SHAPETRIANGLEQUA
#define SHAPETRIANGLEQUA

#include "DataTypes.h"
#include "IntRuleTriangle.h"

class ShapeTriangleQua{
public:
    // Type definition of local integration rule
    typedef IntRuleTriangle LocIntRule;

    // Element dimension
    const static int Dimension = 2;

    // Element number of sides - For Lagrangian elements, this is the same as the number of nodes
    const static int NSides = 6;

    // Number of element nodes
    const static int NElNodes = 6;

    // Number of corner (vertex) nodes
    const static int NCornerNodes = 3;

    const static ElementType ElType = ElementType::ETriangle;

    /// Shapes the shape function value
    /// @param VecDouble Parametric coordinates 
    /// @param VecDouble Shape function values
    static void Shape(VecDouble &xi, VecDouble &phi);
    
    /// Shapes the values of the shape funtion derivatives
    /// @param VecDouble Parametric coordinates 
    /// @param MatrixDouble Shape function derivatives values
    static void ShapeGradient(VecDouble &xi, MatrixDouble &dphi);

    /// Shapes the values of the shape funtion second derivatives    
    /// @param VecDouble Parametric coordinates 
    /// @param std::vector<MatrixDouble> Shape function second derivatives values
    static void ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi);

    /// Gets the element nodal coordinates
    /// @param VecDouble Parametric coordinates 
    static void getCoordinates(MatrixDouble &coord);
   
};

#endif
