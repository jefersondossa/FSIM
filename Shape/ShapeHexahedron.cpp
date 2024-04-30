#include "ShapeHexahedron.h"
const int ShapeHexahedron::Dimension;
const int ShapeHexahedron::Order;
const int ShapeHexahedron::NElNodes;
const int ShapeHexahedron::NCornerNodes;
const ElementType ShapeHexahedron::ElType;

void ShapeHexahedron::Shape(VecDouble &xi, VecDouble &phi) {
    double xsi1 = xi[0];
    double xsi2 = xi[1];
    double xsi3 = xi[2];

    phi[0] = 0.125 * (1. - xsi1)*(1. - xsi2)*(1. - xsi3);
    phi[1] = 0.125 * (1. + xsi1)*(1. - xsi2)*(1. - xsi3);
    phi[2] = 0.125 * (1. + xsi1)*(1. + xsi2)*(1. - xsi3);
    phi[3] = 0.125 * (1. - xsi1)*(1. + xsi2)*(1. - xsi3);
    phi[4] = 0.125 * (1. - xsi1)*(1. - xsi2)*(1. + xsi3);
    phi[5] = 0.125 * (1. + xsi1)*(1. - xsi2)*(1. + xsi3);
    phi[6] = 0.125 * (1. + xsi1)*(1. + xsi2)*(1. + xsi3);
    phi[7] = 0.125 * (1. - xsi1)*(1. + xsi2)*(1. + xsi3);
}

void ShapeHexahedron::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    const double xsi1 = xi[0];
    const double xsi2 = xi[1];
    const double xsi3 = xi[2];
    
    dphi(0,0) = -0.125 * (1. - xsi2) * (1. - xsi3);
    dphi(0,1) = -0.125 * (1. - xsi1) * (1. - xsi3);
    dphi(0,2) = -0.125 * (1. - xsi1) * (1. - xsi2);

    dphi(1,0) =  0.125 * (1. - xsi2) * (1. - xsi3);
    dphi(1,1) = -0.125 * (1. + xsi1) * (1. - xsi3);
    dphi(1,2) = -0.125 * (1. + xsi1) * (1. - xsi2);

    dphi(2,0) =  0.125 * (1. + xsi2) * (1. - xsi3);
    dphi(2,1) =  0.125 * (1. + xsi1) * (1. - xsi3);
    dphi(2,2) = -0.125 * (1. + xsi1) * (1. + xsi2);

    dphi(3,0) = -0.125 * (1. + xsi2) * (1. - xsi3);
    dphi(3,1) =  0.125 * (1. - xsi1) * (1. - xsi3);
    dphi(3,2) = -0.125 * (1. - xsi1) * (1. + xsi2);

    dphi(4,0) = -0.125 * (1. - xsi2) * (1. + xsi3);
    dphi(4,1) = -0.125 * (1. - xsi1) * (1. + xsi3);
    dphi(4,2) =  0.125 * (1. - xsi1) * (1. - xsi2);

    dphi(5,0) =  0.125 * (1. - xsi2) * (1. + xsi3);
    dphi(5,1) = -0.125 * (1. + xsi1) * (1. + xsi3);
    dphi(5,2) =  0.125 * (1. + xsi1) * (1. - xsi2);

    dphi(6,0) =  0.125 * (1. + xsi2) * (1. + xsi3);
    dphi(6,1) =  0.125 * (1. + xsi1) * (1. + xsi3);
    dphi(6,2) =  0.125 * (1. + xsi1) * (1. + xsi2);

    dphi(7,0) = -0.125 * (1. + xsi2) * (1. + xsi3);
    dphi(7,1) =  0.125 * (1. - xsi1) * (1. + xsi3);
    dphi(7,2) =  0.125 * (1. - xsi1) * (1. + xsi2);
}

void ShapeHexahedron::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    return;
}

void ShapeHexahedron::getCoordinates(MatrixDouble &coord) {

    coord(0,0) = -1.0;
    coord(1,0) = -1.0;
    coord(2,0) = -1.0;

    coord(0,1) =  1.0;
    coord(1,1) = -1.0;
    coord(2,1) = -1.0;

    coord(0,2) =  1.0;
    coord(1,2) =  1.0;
    coord(2,2) = -1.0;

    coord(0,3) = -1.0;
    coord(1,3) =  1.0;
    coord(2,3) = -1.0;

    coord(0,4) = -1.0;
    coord(1,4) = -1.0;
    coord(2,4) =  1.0;

    coord(0,5) =  1.0;
    coord(1,5) = -1.0;
    coord(2,5) =  1.0;

    coord(0,6) =  1.0;
    coord(1,6) =  1.0;
    coord(2,6) =  1.0;

    coord(0,7) = -1.0;
    coord(1,7) =  1.0;
    coord(2,7) =  1.0;
}
