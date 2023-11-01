#include "ShapeOneD.h"

const int ShapeOneD::Dimension;

void ShapeOneD::Shape(VecDouble &xi, VecDouble &phi) const {
    PanicButton();
}

void ShapeOneD::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) const {
    PanicButton();

}

void ShapeOneD::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) const {
    PanicButton();

}