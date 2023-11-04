#include "ShapeOneDQua.h"

const int ShapeOneDQua::Dimension;
const int ShapeOneDQua::Order;
const int ShapeOneDQua::NElNodes;

void ShapeOneDQua::Shape(VecDouble &xi, VecDouble &phi) {
    phi[0] = (xi[0] - 1.) * xi[0] / 2.;
    phi[1] = (1. + xi[0]) * xi[0] / 2.;
    phi[2] = (1. - xi[0] * xi[0]); 
}

void ShapeOneDQua::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    dphi(0,0) = xi[0] - 0.5;
    dphi(1,0) = xi[0] + 0.5;
    dphi(2,0) = -2. * xi[0];
}

void ShapeOneDQua::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    std::cout << "Please implement me\n";
    PanicButton();
    return;
}


void ShapeOneDQua::getCoordinates(MatrixDouble &coord) {

    std::cout << "Not implemented\n";
    PanicButton();
}
