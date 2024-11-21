#include "ShapeOneDLin.h"

const int ShapeOneDLin::Dimension;
const int ShapeOneDLin::Order;
const int ShapeOneDLin::NElNodes;
const int ShapeOneDLin::NCornerNodes;
const ElementType ShapeOneDLin::ElType;

void ShapeOneDLin::Shape(VecDouble &xi, VecDouble &phi) {
    phi[0] = (1 - xi[0]) / 2.;
    phi[1] = (1 + xi[0]) / 2.;
}

void ShapeOneDLin::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    dphi(0,0) = -0.5;
    dphi(0,1) =  0.5;
}

void ShapeOneDLin::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    return;
}


void ShapeOneDLin::getCoordinates(MatrixDouble &coord) {
    coord(0,0) = -1.0;
    coord(0,1) =  1.0;
}
