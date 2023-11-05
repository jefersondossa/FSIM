#include "ShapeQuadrilateralLin.h"

const int ShapeQuadrilateralLin::Dimension;
const int ShapeQuadrilateralLin::Order;
const int ShapeQuadrilateralLin::NElNodes;
const int ShapeQuadrilateralLin::NCornerNodes;

void ShapeQuadrilateralLin::Shape(VecDouble &xi, VecDouble &phi) {
    PanicButton();
}

void ShapeQuadrilateralLin::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    PanicButton();

}

void ShapeQuadrilateralLin::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    PanicButton();

}


void ShapeQuadrilateralLin::getCoordinates(MatrixDouble &coord) {

    std::cout << "Not implemented\n";
    PanicButton();
}