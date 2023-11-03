#include "ShapeOneD.h"

const int ShapeOneD::Dimension;
const int ShapeOneD::Order;
const int ShapeOneD::NElNodes;

void ShapeOneD::Shape(VecDouble &xi, VecDouble &phi) {
    phi[0] = (1 - xi[0]) / 2.;
    phi[1] = (1 + xi[0]) / 2.;
}

void ShapeOneD::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    dphi(0,0) = -0.5;
    dphi(1,0) =  0.5;
}

void ShapeOneD::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    return;
}


void ShapeOneD::getCoordinates(MatrixDouble &coord) {

    std::cout << "Not implemented\n";
    PanicButton();
}
