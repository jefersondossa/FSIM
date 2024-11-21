#include "ShapeTriangleLin.h"

const int ShapeTriangleLin::Dimension;
const int ShapeTriangleLin::Order;
const int ShapeTriangleLin::NElNodes;
const int ShapeTriangleLin::NCornerNodes;
const ElementType ShapeTriangleLin::ElType;

void ShapeTriangleLin::Shape(VecDouble &xi, VecDouble &phi) {

    double xsi1 = xi[0];
    double xsi2 = xi[1];
    double xsi3 = 1. - xsi1 - xsi2;

    phi[0] = xsi3;
    phi[1] = xsi1;
    phi[2] = xsi2;
    // element conectivity
    //     2
    //     01

}

void ShapeTriangleLin::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    const double xsi1 = xi[0];
    const double xsi2 = xi[1];
    const double xsi3 = 1. - xsi1 - xsi2;

    dphi(0,0) = -1.;
    dphi(1,0) = -1.;
    dphi(0,1) = 1.;
    dphi(1,1) = 0.;
    dphi(0,2) = 0.;
    dphi(1,2) = 1.;
    // element conectivity
    //     2
    //     01
}

void ShapeTriangleLin::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    
    return;
}


void ShapeTriangleLin::getCoordinates(MatrixDouble &coord) {

    coord(0,0) = 0.0;
    coord(1,0) = 0.0;

    coord(0,1) = 1.0;
    coord(1,1) = 0.0;

    coord(0,2) = 0.0;
    coord(1,2) = 1.0;

}
