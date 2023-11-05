#ifndef SHAPETETRAHEDRONCUB
#define SHAPETETRAHEDRONCUB

#include "DataTypes.h"
#include "IntRuleTetrahedron.h"

class ShapeTetrahedronCub{
public:
    // Type definition of integration rule for quadrilateral elements
    typedef IntRuleTetrahedron LocIntRule;

    const static int Dimension = 3;
    const static int Order = 3;
    const static int NElNodes = 20;
    const static int NCornerNodes = 4;
    
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
