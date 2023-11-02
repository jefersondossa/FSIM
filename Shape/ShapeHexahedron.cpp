#include "ShapeHexahedron.h"
const int ShapeHexahedron::Dimension;
const int ShapeHexahedron::Order;
const int ShapeHexahedron::NElNodes;


void ShapeHexahedron::Shape(VecDouble &xi, VecDouble &phi) {
    PanicButton();
}

void ShapeHexahedron::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    PanicButton();

}

void ShapeHexahedron::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    PanicButton();

}

void ShapeHexahedron::getCoordinates(MatrixDouble &coord) {

    std::cout << "Not implemented\n";
    PanicButton();
}
