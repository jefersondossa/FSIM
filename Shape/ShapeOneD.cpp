#include "ShapeOneD.h"

const int ShapeOneD::Dimension;
const int ShapeOneD::Order;
const int ShapeOneD::NElNodes;

void ShapeOneD::Shape(VecDouble &xi, VecDouble &phi) {
    PanicButton();
}

void ShapeOneD::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    PanicButton();

}

void ShapeOneD::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    PanicButton();

}


void ShapeOneD::getCoordinates(MatrixDouble &coord) {

    std::cout << "Not implemented\n";
    PanicButton();
}
