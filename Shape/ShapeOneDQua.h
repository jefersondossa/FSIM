#ifndef ShapeOneDQua_H
#define ShapeOneDQua_H

#include "DataTypes.h"
#include "IntRule1d.h"

class ShapeOneDQua{
public:
    // Type definition of local integration rule
    typedef IntRule1d LocIntRule;

    // Element dimension
    const static int Dimension = 1;

    // Element number of sides - For Lagrangian elements, this is the same as the number of nodes
    const static int NSides = 3;

    // Number of element nodes
    const static int NShape = 3;

    // Number of corner (vertex) nodes
    const static int NCornerNodes = 2;

    const static ElementType ElType = ElementType::EOneD;

    static int NShapeFunctions(int side, int order){return NShape;};
    static int NShapeFunctions(int order){return NShape;};
    
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
