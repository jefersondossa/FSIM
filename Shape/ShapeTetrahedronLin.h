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

    // Element number of sides - For Lagrangian elements, this is the same as the number of nodes
    const static int NSides = 4;

    // Number of element nodes
    const static int NShape = 4;

    // Number of corner (vertex) nodes
    const static int NCornerNodes = 4;
    const static int NEdges = 6;
    const static int NFaces = 4;
    const static int NVolumes = 1;

    const static ElementType ElType = ElementType::ETetrahedron;

    static int NShapeFunctions(int side, int order){return 1;};
    static int NShapeFunctions(int order){return NShape;};

    /// Shapes the shape function value
    /// @param VecDouble Parametric coordinates 
    /// @param VecDouble Shape function values
    /// @param MatrixDouble Shape function gradient values
    static void Shape(VecDouble &xi, VecDouble &phi, MatrixDouble &dphi, VecInt orders);

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
