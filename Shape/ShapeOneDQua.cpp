#include "ShapeOneDQua.h"

const int ShapeOneDQua::Dimension;
const int ShapeOneDQua::NSides;
const int ShapeOneDQua::NShape;
const int ShapeOneDQua::NCornerNodes;
const ElementType ShapeOneDQua::ElType;

void ShapeOneDQua::Shape(VecDouble &xi, VecDouble &phi) {
    phi[0] = (xi[0] - 1.) * xi[0] / 2.;
    phi[1] = (1. + xi[0]) * xi[0] / 2.;
    phi[2] = (1. - xi[0] * xi[0]); 
}

void ShapeOneDQua::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    dphi(0,0) = xi[0] - 0.5;
    dphi(0,1) = xi[0] + 0.5;
    dphi(0,2) = -2. * xi[0];
}

void ShapeOneDQua::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    ddphi[0](0,0) = 1.;
    ddphi[1](0,0) = 1.;
    ddphi[2](0,0) =-2.;
}


void ShapeOneDQua::getCoordinates(MatrixDouble &coord) {
    coord(0,0) = -1.0;
    coord(0,1) =  1.0;
    coord(0,2) =  0.0;
}
