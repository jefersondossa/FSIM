#include "ShapeQuadrilateralLin.h"

const int ShapeQuadrilateralLin::Dimension;
const int ShapeQuadrilateralLin::Order;
const int ShapeQuadrilateralLin::NElNodes;
const int ShapeQuadrilateralLin::NCornerNodes;

void ShapeQuadrilateralLin::Shape(VecDouble &xi, VecDouble &phi) {
    phi[0] = (1 - xi[0])*(1 - xi[1]) / 4.;
    phi[1] = (1 + xi[0])*(1 - xi[1]) / 4.;
    phi[2] = (1 + xi[0])*(1 + xi[1]) / 4.;
    phi[3] = (1 - xi[0])*(1 + xi[1]) / 4.;
}

void ShapeQuadrilateralLin::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    dphi(0,0) = -0.25 * (1. - xi[1]);
    dphi(0,1) = -0.25 * (1. - xi[0]);

    dphi(1,0) =  0.25 * (1. - xi[1]);
    dphi(1,1) = -0.25 * (1. + xi[0]);

    dphi(2,0) =  0.25 * (1. + xi[1]);
    dphi(2,1) =  0.25 * (1. + xi[0]);

    dphi(3,0) = -0.25 * (1. + xi[1]);
    dphi(3,1) =  0.25 * (1. - xi[0]);
}

void ShapeQuadrilateralLin::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    return;
}


void ShapeQuadrilateralLin::getCoordinates(MatrixDouble &coord) {
    coord(0,0) = -1.0;
    coord(1,0) = -1.0;

    coord(0,1) = -1.0;
    coord(1,1) =  1.0;

    coord(0,2) =  1.0;
    coord(1,2) =  1.0;

    coord(0,3) =  1.0;
    coord(1,3) = -1.0;
    
}