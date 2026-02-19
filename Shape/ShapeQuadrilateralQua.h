#ifndef SHAPEQUADRILATERALQUA
#define SHAPEQUADRILATERALQUA

#include "DataTypes.h"
#include "IntRuleQuad.h"

class ShapeQuadrilateralQua{
public:
    // Type definition of local integration rule
    typedef IntRuleQuad LocIntRule;

    // Element dimension
    const static int Dimension = 2;

    // Element number of sides - For Lagrangian elements, this is the same as the number of nodes
    const static int NSides = 9;

    // Number of element nodes
    const static int NShape = 9;

    // Number of corner (vertex) nodes
    const static int NCornerNodes = 4;
    const static int NEdges = 4;
    const static int NFaces = 1;
    const static int NVolumes = 0;

    const static ElementType ElType = ElementType::EQuadrilateral;
    
    static int NShapeFunctions(int side, int order){return 1;};
    static int NShapeFunctions(int order){return NShape;};
    
    /// Shapes the shape function value
    /// @param VecDouble Parametric coordinates 
    /// @param VecDouble Shape function values
    static void Shape(VecDouble &xi, VecDouble &phi, MatrixDouble &dphi, int order = 1);

    /// Shapes the values of the shape funtion second derivatives    
    /// @param VecDouble Parametric coordinates 
    /// @param std::vector<MatrixDouble> Shape function second derivatives values
    static void ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi);

    /// Gets the element nodal coordinates
    /// @param VecDouble Parametric coordinates 
    static void getCoordinates(MatrixDouble &coord);
   
    // Number of nodes associated with a side
    static int NSideNodes(int side){return 1;};

    static int SideNodeLocIndex(int side, int node){
        return node;
    };
};

#endif
