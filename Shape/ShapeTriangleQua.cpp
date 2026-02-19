#include "ShapeTriangleQua.h"

const int ShapeTriangleQua::Dimension;
const int ShapeTriangleQua::NSides;
const int ShapeTriangleQua::NShape;
const int ShapeTriangleQua::NCornerNodes;
const int ShapeTriangleQua::NEdges;
const int ShapeTriangleQua::NFaces;
const int ShapeTriangleQua::NVolumes;
const ElementType ShapeTriangleQua::ElType;

void ShapeTriangleQua::Shape(VecDouble &xi, VecDouble &phi, MatrixDouble &dphi, int order) {

    double xsi1 = xi[0];
    double xsi2 = xi[1];
    double xsi3 = 1. - xsi1 - xsi2;
    
    phi[0] = xsi3 * (2.0 * xsi3 - 1.0);
    phi[1] = xsi1 * (2.0 * xsi1 - 1.0);
    phi[2] = xsi2 * (2.0 * xsi2 - 1.0);
    phi[3] = 4.0 * xsi3 * xsi1;
    phi[4] = 4.0 * xsi1 * xsi2;
    phi[5] = 4.0 * xsi2 * xsi3;

    dphi(0,1) = 4. * xsi1 - 1.;
    dphi(1,1) = 0.;
    dphi(0,2) = 0.;
    dphi(1,2) = 4. * xsi2 - 1.;
    dphi(0,0) = -4. * xsi3 + 1.;
    dphi(1,0) = -4. * xsi3 + 1.;
    dphi(0,4) = 4. * xsi2;
    dphi(1,4) = 4. * xsi1;
    dphi(0,5) = -4. * xsi2;
    dphi(1,5) = 4. * (xsi3 - xsi2);
    dphi(0,3) = 4. * (xsi3 - xsi1);
    dphi(1,3) = -4. * xsi1;
    // element conectivity
    //     2
    //     54
    //     031

}

void ShapeTriangleQua::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    const double xsi1 = xi[0];
    const double xsi2 = xi[1];
    
    ddphi[0](0,0) = 4.;
    ddphi[0](0,1) = 4.;
    ddphi[0](1,0) = 4.;
    ddphi[0](1,1) = 4.;
    ddphi[1](0,0) = 4.;
    ddphi[1](0,1) = 0.;
    ddphi[1](1,0) = 0.;
    ddphi[1](1,1) = 0.;
    ddphi[2](0,0) = 0.;
    ddphi[2](0,1) = 0.;
    ddphi[2](1,0) = 0.;
    ddphi[2](1,1) = 4.;
    ddphi[3](0,0) = -8.;
    ddphi[3](0,1) = -4.;
    ddphi[3](1,0) = -4.;
    ddphi[3](1,1) = 0.;
    ddphi[4](0,0) = 0.;
    ddphi[4](0,1) = 4.;
    ddphi[4](1,0) = 4.;
    ddphi[4](1,1) = 0.;
    ddphi[5](0,0) = 0.;
    ddphi[5](0,1) = -4.;
    ddphi[5](1,0) = -4.;
    ddphi[5](1,1) = -8.;
    // element conectivity
    //     2
    //     54
    //     031
}


void ShapeTriangleQua::getCoordinates(MatrixDouble &coord) {

    coord(0,0) = 0.0;
    coord(1,0) = 0.0;

    coord(0,1) = 1.0;
    coord(1,1) = 0.0;

    coord(0,2) = 0.0;
    coord(1,2) = 1.0;

    coord(0,3) = 0.5;
    coord(1,3) = 0.0;

    coord(0,4) = 0.5;
    coord(1,4) = 0.5;

    coord(0,5) = 0.0;
    coord(1,5) = 0.5;

}
