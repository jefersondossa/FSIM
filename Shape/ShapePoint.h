#ifndef SHAPEPOINT
#define SHAPEPOINT

#include "DataTypes.h"
#include "IntRule0d.h"

class ShapePoint{
public:
    // Type definition for integration rule 1D
    typedef IntRule0d LocIntRule;
    
    const static int Dimension = 0; 
    const static int Order = 1;
    const static int NElNodes = 1;

    /// Shapes the shape function value
    /// @param double* Non-dimensional coordinates 
    /// @param double* Shape function values
    static void Shape(VecDouble &xi, VecDouble &phi){
        phi[0] = 1.;
    };
    
    /// Shapes the values of the shape funtion derivatives
    /// @param double* Non-dimensional coordinates 
    /// @param double** Shape function derivatives values
    static void ShapeGradient(VecDouble &xi, MatrixDouble &dphi){
        dphi.setZero();
    };   

    /// Shapes the values of the shape funtion second derivatives    
    /// @param double** Shape function second derivatives values
    static void ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi){
        return;
    };     

    /// Gets the element nodal coordinates
    /// @param double** Non-dimensional coordinate
    static void getCoordinates(MatrixDouble &coord){
        coord(0,0) = 0.;
        coord(1,0) = 0.;
        coord(2,0) = 0.;
    };


    
};

#endif
