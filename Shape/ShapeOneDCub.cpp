#include "ShapeOneDCub.h"

const int ShapeOneDCub::Dimension;
const int ShapeOneDCub::Order;
const int ShapeOneDCub::NElNodes;

void ShapeOneDCub::Shape(VecDouble &xi, VecDouble &phi) {
    phi[0] = (-1. + xi[0] + 9.*xi[0]*xi[0] - 9.*xi[0]*xi[0]*xi[0]) / 16.;
    phi[1] = (1. + xi[0])*(-1. + 3.*xi[0])*(1. + 3.*xi[0]) / 16.;
    phi[2] = (1. + xi[0])*(1. - 4.*xi[0] + 3.*xi[0]*xi[0]) * 9. / 16.;
    phi[3] = (-1. + xi[0])*(1. + xi[0])*(1.+3.*xi[0]) * - 9. / 16.;
}

void ShapeOneDCub::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    dphi(0,0) = (1. + 18. * xi[0] - 27. * xi[0]*xi[0])/16.;
    dphi(1,0) = (-1. + 18. * xi[0] + 27. * xi[0]*xi[0])/16.;
    dphi(2,0) = (-3. - 2. * xi[0] + 9. * xi[0]*xi[0])*9./16.;
    dphi(3,0) = -(-3. + 2. * xi[0] + 9. * xi[0]*xi[0])*9./16.;
}

void ShapeOneDCub::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    return;
}


void ShapeOneDCub::getCoordinates(MatrixDouble &coord) {

    std::cout << "Not implemented\n";
    PanicButton();
}
