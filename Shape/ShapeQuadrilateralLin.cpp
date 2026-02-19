#include "ShapeQuadrilateralLin.h"

const int ShapeQuadrilateralLin::Dimension;
const int ShapeQuadrilateralLin::NSides;
const int ShapeQuadrilateralLin::NShape;
const int ShapeQuadrilateralLin::NCornerNodes;
const int ShapeQuadrilateralLin::NEdges;
const int ShapeQuadrilateralLin::NFaces;
const int ShapeQuadrilateralLin::NVolumes;
const ElementType ShapeQuadrilateralLin::ElType;

void ShapeQuadrilateralLin::Shape(VecDouble &xi, VecDouble &phi, MatrixDouble &dphi, VecInt orders) {
    phi[0] = 0.25 * (1 - xi[0])*(1 - xi[1]);
    phi[1] = 0.25 * (1 + xi[0])*(1 - xi[1]);
    phi[2] = 0.25 * (1 + xi[0])*(1 + xi[1]);
    phi[3] = 0.25 * (1 - xi[0])*(1 + xi[1]);

    dphi(0,0) = -0.25 * (1. - xi[1]);
    dphi(1,0) = -0.25 * (1. - xi[0]);
    dphi(0,1) =  0.25 * (1. - xi[1]);
    dphi(1,1) = -0.25 * (1. + xi[0]);
    dphi(0,2) =  0.25 * (1. + xi[1]);
    dphi(1,2) =  0.25 * (1. + xi[0]);
    dphi(0,3) = -0.25 * (1. + xi[1]);
    dphi(1,3) =  0.25 * (1. - xi[0]);
}

void ShapeQuadrilateralLin::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    const double xsi1 = xi[0];
    const double xsi2 = xi[1];
    const double xsi3 = xi[2];

    ddphi[0](0,0) =  0.;
    ddphi[0](0,1) =  0.25;
    ddphi[0](1,0) =  0.25;
    ddphi[0](1,1) =  0.;
    
    ddphi[1](0,0) =  0.;
    ddphi[1](0,1) = -0.25;
    ddphi[1](1,0) = -0.25;
    ddphi[1](1,1) =  0.;

    ddphi[2](0,0) =  0.;
    ddphi[2](0,1) =  0.25;
    ddphi[2](1,0) =  0.25;
    ddphi[2](1,1) =  0.;

    ddphi[3](0,0) =  0.;
    ddphi[3](0,1) = -0.25;
    ddphi[3](1,0) = -0.25;
    ddphi[3](1,1) =  0.;
    
    return;
}


void ShapeQuadrilateralLin::getCoordinates(MatrixDouble &coord) {
    coord(0,0) = -1.0;
    coord(1,0) = -1.0;

    coord(0,1) =  1.0;
    coord(1,1) = -1.0;

    coord(0,2) =  1.0;
    coord(1,2) =  1.0;

    coord(0,3) = -1.0;
    coord(1,3) =  1.0;
    
}