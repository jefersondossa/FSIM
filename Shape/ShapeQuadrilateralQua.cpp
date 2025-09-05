#include "ShapeQuadrilateralQua.h"

const int ShapeQuadrilateralQua::Dimension;
const int ShapeQuadrilateralQua::NSides;
const int ShapeQuadrilateralQua::NElNodes;
const int ShapeQuadrilateralQua::NCornerNodes;
const ElementType ShapeQuadrilateralQua::ElType;

void ShapeQuadrilateralQua::Shape(VecDouble &xi, VecDouble &phi) {
    const double xsi1 = xi[0];
    const double xsi2 = xi[1];
    
    phi[0] = 0.25 * (xsi1 - 1.) * (xsi2 - 1.) * xsi1 * xsi2;
    phi[1] = 0.25 * (xsi1 + 1.) * (xsi2 - 1.) * xsi1 * xsi2;
    phi[2] = 0.25 * (xsi1 + 1.) * (xsi2 + 1.) * xsi1 * xsi2;
    phi[3] = 0.25 * (xsi1 - 1.) * (xsi2 + 1.) * xsi1 * xsi2;
    phi[4] =  0.5 * (1. - xsi1*xsi1) * (xsi2 - 1.) * xsi2;
    phi[5] =  0.5 * xsi1 * (xsi1 + 1.) * (1. - xsi2*xsi2);
    phi[6] =  0.5 * (1. - xsi1*xsi1) * (xsi2 + 1.) * xsi2;
    phi[7] =  0.5 * xsi1 * (xsi1 - 1.) * (1. - xsi2*xsi2);
    phi[8] =  (1. - xsi1*xsi1) * (1. - xsi2*xsi2);
}

void ShapeQuadrilateralQua::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    const double xsi1 = xi[0];
    const double xsi2 = xi[1];

    dphi(0,0) = 0.5 * (-0.5 + 1. * xsi1) * (-1. + xsi2) * xsi2;
    dphi(1,0) = 0.5 * (-1. + xsi1) * xsi1 * (-0.5 + 1. * xsi2);
    dphi(0,1) = 0.5 * (0.5 + 1. * xsi1) * (-1. + xsi2) * xsi2;
    dphi(1,1) = 0.5 * xsi1 * (1. + xsi1) * (-0.5 + xsi2);
    dphi(0,2) = 0.5 * (0.5 + 1. * xsi1) * xsi2 * (1. + xsi2);
    dphi(1,2) = 0.5 * xsi1 * (1. + xsi1) * (0.5 + xsi2);
    dphi(0,3) = 0.5 * (-0.5 + xsi1) * xsi2 * (1. + xsi2);
    dphi(1,3) = 0.5 * (-1. + xsi1) * xsi1 * (0.5 + xsi2);
    dphi(0,4) = -xsi1 * (-1. + xsi2) * xsi2;
    dphi(1,4) = -0.5 + xsi1 * xsi1 * (0.5 - xsi2) + xsi2;
    dphi(0,5) = 0.5 - 0.5 * xsi2 * xsi2 + xsi1 * (1. - xsi2 * xsi2);
    dphi(1,5) = -xsi1 * (1. + xsi1) * xsi2;
    dphi(0,6) = -xsi1 * xsi2 * (1. + xsi2);
    dphi(1,6) = 0.5 + xsi1 * xsi1 * (-0.5 - xsi2) + xsi2;
    dphi(0,7) = -0.5 + 0.5 * xsi2 * xsi2 + xsi1 * (1. - xsi2 * xsi2);
    dphi(1,7) = -(-1. + xsi1) * xsi1 * xsi2;
    dphi(0,8) = 2. * xsi1 * (-1. + xsi2 * xsi2);
    dphi(1,8) = 2. * (-1. + xsi1 * xsi1) * xsi2;
}

void ShapeQuadrilateralQua::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    std::cout << "Please implement me " << std::endl;
    PanicButton();
    return;
}


void ShapeQuadrilateralQua::getCoordinates(MatrixDouble &coord) {
    coord(0,0) = -1.0;
    coord(1,0) = -1.0;

    coord(0,1) =  1.0;
    coord(1,1) = -1.0;

    coord(0,2) =  1.0;
    coord(1,2) =  1.0;

    coord(0,3) = -1.0;
    coord(1,3) =  1.0;

    coord(0,4) =  0.0;
    coord(1,4) = -1.0;

    coord(0,5) =  1.0;
    coord(1,5) =  0.0;

    coord(0,6) =  0.0;
    coord(1,6) =  1.0;

    coord(0,7) = -1.0;
    coord(1,7) =  0.0;

    coord(0,8) =  0.0;
    coord(1,8) =  0.0;
    
}