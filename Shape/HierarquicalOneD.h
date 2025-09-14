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

    constexpr static int EdgePermutation[2][2] = {{0,1},{1,0}};

    static VecInt GetEdgePermutation(int side){
        if(side < 0 || side > 1) PanicButton();
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
    static int NShapeFunctions(int order){return order+1;};
    
    // Returns the total number of shape functions
    static int NShapeFunctions(VecInt &orders);

    /// Shapes the shape function value
    /// @param VecDouble Parametric coordinates 
    /// @param VecDouble Shape function values
    static void Shape(VecDouble &xi, VecDouble &phi, int order);
    
    /// Shapes the values of the shape funtion derivatives
    /// @param VecDouble Parametric coordinates 
    /// @param MatrixDouble Shape function derivatives values
    static void ShapeGradient(VecDouble &xi, MatrixDouble &dphi, int order);

    /// Shapes the values of the shape funtion second derivatives    
    /// @param VecDouble Parametric coordinates 
    /// @param std::vector<MatrixDouble> Shape function second derivatives values
    static void ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi);  



    static void Chebyshev(double x, int num,MatrixDouble &phi,MatrixDouble &dphi);
};

#endif
