#ifndef SHAPEPOINT
#define SHAPEPOINT

#include "ShapeFunction.h"

class ShapePoint : public ShapeFunction{
public:
    const static int Dimension = 0; 

    /// Shapes the shape function value
    /// @param double* Non-dimensional coordinates 
    /// @param double* Shape function values
    void Shape(VecDouble &xi, VecDouble &phi) const override{
        phi[0] = 1.;
    };
    
    /// Shapes the values of the shape funtion derivatives
    /// @param double* Non-dimensional coordinates 
    /// @param double** Shape function derivatives values
    void ShapeGradient(VecDouble &xi, MatrixDouble &dphi) const override{
        dphi.setZero();
    };   

    /// Shapes the values of the shape funtion second derivatives    
    /// @param double** Shape function second derivatives values
    void ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) const override{
        return;
    };     

    /// Gets the element nodal coordinates
    /// @param double** Non-dimensional coordinate
    void getCoordinates(MatrixDouble &coord) const override{
        coord(0,0) = 0.;
        coord(1,0) = 0.;
        coord(2,0) = 0.;
    };


    
};

#endif
