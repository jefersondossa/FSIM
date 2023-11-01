#include "ShapeQuadrilateral.h"
const int ShapeQuadrilateral::Dimension;

void ShapeQuadrilateral::Shape(VecDouble &xi, VecDouble &phi) const {
    PanicButton();
}

void ShapeQuadrilateral::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) const {
    PanicButton();

}

void ShapeQuadrilateral::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) const {
    PanicButton();

}