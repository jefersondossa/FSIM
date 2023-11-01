#ifndef SHAPEHEXAHEDRON
#define SHAPEHEXAHEDRON

#include "ShapeFunction.h"

class ShapeHexahedron : public ShapeFunction{
public:
    const static int Dimension = 3;

    /// Shapes the shape function value
    /// @param double* Non-dimensional coordinates 
    /// @param double* Shape function values
    void Shape(VecDouble &xi, VecDouble &phi) const override;
    
    /// Shapes the values of the shape funtion derivatives
    /// @param double* Non-dimensional coordinates 
    /// @param double** Shape function derivatives values
    void ShapeGradient(VecDouble &xi, MatrixDouble &dphi) const override;

    /// Shapes the values of the shape funtion second derivatives    
    /// @param double** Shape function second derivatives values
    void ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) const override;

    /// Gets the element nodal coordinates
    /// @param double** Non-dimensional coordinate
    void getCoordinates(MatrixDouble &coord) const override;
};

#endif
