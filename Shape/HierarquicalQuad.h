#ifndef HierarquicalQuad_H
#define HierarquicalQuad_H

#include "DataTypes.h"
#include "IntRuleQuad.h"

class HierarquicalQuad{
public:
    // Type definition of local integration rule
    typedef IntRuleQuad LocIntRule;

    // Element dimension
    const static int Dimension = 2;

    // Number of corner (vertex) nodes
    const static int NCornerNodes = 4;

    const static int NSides = 9;

    const static ElementType ElType = ElementType::EQuadrilateral;

    constexpr static int EdgePermutation[8][2] = {{0,1},{1,0},{1,2},{2,1},{2,3},{3,2},{3,0},{0,3}};

    static VecInt GetEdgePermutation(int side){
        if(side < 0 || side > 7) PanicButton();
        VecInt perm(2);
        perm[0] = EdgePermutation[side][0];
        perm[1] = EdgePermutation[side][1];
        return perm;
    };

    // Number of nodes associated with a side
    static int NSideNodes(int side);
    
    // Local node index of a node associated with a side
    static int SideNodeLocIndex(int side, int node);

    // Returns the number of shape functions associated with a side
    static int NShapeFunctions(int side, int order);
    static int NShapeFunctions(int order){return (order+1)*(order+1);};
    
    // Returns the total number of shape functions
    static int NShapeFunctions(VecInt &orders);

    /// Shapes the shape function value
    /// @param VecDouble Parametric coordinates 
    /// @param VecDouble Shape function values
    static void Shape(VecDouble &xi, VecDouble &phi, MatrixDouble &dphi, VecInt orders);
    
    /// Shapes the values of the shape funtion second derivatives    
    /// @param VecDouble Parametric coordinates 
    /// @param std::vector<MatrixDouble> Shape function second derivatives values
    static void ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi);  

    static void Chebyshev(REAL x, int num,MatrixDouble &phi,MatrixDouble &dphi);
    static void Legendre(REAL x, int num,MatrixDouble &phi,MatrixDouble &dphi);
};

#endif
