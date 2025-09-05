#ifndef HierarquicalOneD_H
#define HierarquicalOneD_H

#include "DataTypes.h"
#include "IntRule1d.h"

class HierarquicalOneD{
public:
    // Type definition of local integration rule
    typedef IntRule1d LocIntRule;

    // Element dimension
    const static int Dimension = 1;

    // Number of corner (vertex) nodes
    const static int NCornerNodes = 2;

    const static int NSides = 3;

    const static ElementType ElType = ElementType::EOneD;
    
    // Number of nodes associated with a side
    static int NSideNodes(int side);
    
    // Local node index of a node associated with a side
    static int SideNodeLocIndex(int side, int node);

    // Returns the number of shape functions associated with a side
    static int NShapeFunctions(int side, int order);
    
    // Returns the total number of shape functions
    static int NShapeFunctions(VecInt &orders);

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
};

#endif
