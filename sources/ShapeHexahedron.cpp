#include "ShapeHexahedron.h"
const int ShapeHexahedron::Dimension;

void ShapeHexahedron::Shape(VecDouble &xi, VecDouble &phi) const {
    PanicButton();
}

void ShapeHexahedron::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) const {
    PanicButton();

}

void ShapeHexahedron::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) const {
    PanicButton();

}