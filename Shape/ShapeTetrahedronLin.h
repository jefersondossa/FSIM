#ifndef SHAPETETRAHEDRONLIN
#define SHAPETETRAHEDRONLIN

#include "DataTypes.h"
#include "IntRuleTetrahedron.h"

class ShapeTetrahedronLin{
public:
    // Type definition of local integration rule
    typedef IntRuleTetrahedron LocIntRule;
    
    // Element dimension
    const static int Dimension = 3;

    // Element polynomial order
    const static int Order = 1;

    // Number of element nodes
    const static int NElNodes = 4;

    // Number of corner (vertex) nodes
    const static int NCornerNodes = 4;

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
