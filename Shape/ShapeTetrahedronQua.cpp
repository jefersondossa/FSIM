#include "ShapeTetrahedronQua.h"

const int ShapeTetrahedronQua::Dimension;
const int ShapeTetrahedronQua::NSides;
const int ShapeTetrahedronQua::NShape;
const int ShapeTetrahedronQua::NCornerNodes;
const int ShapeTetrahedronQua::NEdges;
const int ShapeTetrahedronQua::NFaces;
const int ShapeTetrahedronQua::NVolumes;
const ElementType ShapeTetrahedronQua::ElType;

void ShapeTetrahedronQua::Shape(VecDouble &xi, VecDouble &phi, MatrixDouble &dphi, VecInt orders) {

    REAL xsi1 = xi[0];
    REAL xsi2 = xi[1];
    REAL xsi3 = xi[2];
    REAL xsi4 = 1.0 - xsi1 - xsi2 - xsi3;

    phi[3] = 2.0 * (xsi1 - 0.50) * xsi1;
    phi[1] = 2.0 * (xsi2 - 0.50) * xsi2;
    phi[2] = 2.0 * (xsi3 - 0.50) * xsi3;
    phi[0] = (2.0 - 2.0 * xsi3 - 2.0 * xsi2 - 2.0 * xsi1 - 1.0) * (1.0 - xsi1 - xsi2 - xsi3);
    phi[8] = 4.0 * xsi1 * xsi2;
    phi[9] = 4.0 * xsi1 * xsi3;
    phi[7] = 4.0 * xsi1 * (1.0 - xsi1 - xsi2 - xsi3);
    phi[6] = 4.0 * xsi3 * (1.0 - xsi1 - xsi2 - xsi3);
    phi[4] = 4.0 * xsi2 * (1.0 - xsi1 - xsi2 - xsi3);
    phi[5] = 4.0 * xsi2 * xsi3;

    dphi(0,3) = 4. * xsi1 - 1.;
    dphi(1,3) = 0.;
    dphi(2,3) = 0.;

    dphi(0,1) = 0.;
    dphi(1,1) = 4. * xsi2 - 1.;
    dphi(2,1) = 0.;

    dphi(0,2) = 0.;
    dphi(1,2) = 0.;
    dphi(2,2) = 4. * xsi3 - 1.;

    dphi(0,0) = 4. * (xsi1 + xsi2 + xsi3) - 3.;
    dphi(1,0) = 4. * (xsi1 + xsi2 + xsi3) - 3.;
    dphi(2,0) = 4. * (xsi1 + xsi2 + xsi3) - 3.;

    dphi(0,8) = 4. * xsi2;
    dphi(1,8) = 4. * xsi1;
    dphi(2,8) = 0.;

    dphi(0,9) = 4. * xsi3;
    dphi(1,9) = 0.;
    dphi(2,9) = 4. * xsi1;

    dphi(0,7) = 4. * (1. - 2. * xsi1 - xsi2 - xsi3);
    dphi(1,7) = -4. * xsi1;
    dphi(2,7) = -4. * xsi1;

    dphi(0,5) = 0.;
    dphi(1,5) = 4. * xsi3;
    dphi(2,5) = 4. * xsi2;

    dphi(0,6) = -4. * xsi3;
    dphi(1,6) = -4. * xsi3;
    dphi(2,6) = 4. * (1. - 2. * xsi3 - xsi2 - xsi1);

    dphi(0,4) = -4. * xsi2;
    dphi(1,4) = 4. * (1. - 2. * xsi2 - xsi1 - xsi3);
    dphi(2,4) = -4. * xsi2;   

    // element conectivity
    //layer 1
    //     3
    //     8 9
    //     2 7 1
    //layer 2
    //     6
    //     54
    //layer 3
    //     0
      
}

void ShapeTetrahedronQua::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    const REAL xsi1 = xi[0];
    const REAL xsi2 = xi[1];
    const REAL xsi3 = xi[2];

    ddphi[0](0,0) = 4.;
    ddphi[0](0,1) = 0.;
    ddphi[0](0,2) = 0.;
    ddphi[0](1,0) = 0.;
    ddphi[0](1,1) = 0.;
    ddphi[0](1,2) = 0.;
    ddphi[0](2,0) = 0.;
    ddphi[0](2,1) = 0.;
    ddphi[0](2,2) = 0.;

    ddphi[1](0,0) = 0.;
    ddphi[1](0,1) = 0.;
    ddphi[1](0,2) = 0.;
    ddphi[1](1,0) = 0.;
    ddphi[1](1,1) = 4.;
    ddphi[1](1,2) = 0.;
    ddphi[1](2,0) = 0.;
    ddphi[1](2,1) = 0.;
    ddphi[1](2,2) = 0.;

    ddphi[2](0,0) = 0.;
    ddphi[2](0,1) = 0.;
    ddphi[2](0,2) = 0.;
    ddphi[2](1,0) = 0.;
    ddphi[2](1,1) = 0.;
    ddphi[2](1,2) = 0.;
    ddphi[2](2,0) = 0.;
    ddphi[2](2,1) = 0.;
    ddphi[2](2,2) = 4.;

    ddphi[3](0,0) = 4.;
    ddphi[3](0,1) = 4.;
    ddphi[3](0,2) = 4.;
    ddphi[3](1,0) = 4.;
    ddphi[3](1,1) = 4.;
    ddphi[3](1,2) = 4.;
    ddphi[3](2,0) = 4.;
    ddphi[3](2,1) = 4.;
    ddphi[3](2,2) = 4.;

    ddphi[4](0,0) = 0.;
    ddphi[4](0,1) = 4.;
    ddphi[4](0,2) = 0.;
    ddphi[4](1,0) = 4.;
    ddphi[4](1,1) = 0.;
    ddphi[4](1,2) = 0.;
    ddphi[4](2,0) = 0.;
    ddphi[4](2,1) = 0.;
    ddphi[4](2,2) = 0.;

    ddphi[5](0,0) = 0.;
    ddphi[5](0,1) = 0.;
    ddphi[5](0,2) = 4.;
    ddphi[5](1,0) = 0.;
    ddphi[5](1,1) = 0.;
    ddphi[5](1,2) = 0.;
    ddphi[5](2,0) = 4.;
    ddphi[5](2,1) = 0.;
    ddphi[5](2,2) = 0.;

    ddphi[6](0,0) = -8.;
    ddphi[6](0,1) = -4.;
    ddphi[6](0,2) = -4.;
    ddphi[6](1,0) = -4.;
    ddphi[6](1,1) = 0.;
    ddphi[6](1,2) = 0.;
    ddphi[6](2,0) = -4.;
    ddphi[6](2,1) = 0.;
    ddphi[6](2,2) = 0.;

    ddphi[7](0,0) = 0.;
    ddphi[7](0,1) = 0.;
    ddphi[7](0,2) = 0.;
    ddphi[7](1,0) = 0.;
    ddphi[7](1,1) = 0.;
    ddphi[7](1,2) = 4.;
    ddphi[7](2,0) = 0.;
    ddphi[7](2,1) = 4.;
    ddphi[7](2,2) = 0.;

    ddphi[8](0,0) = 0.;
    ddphi[8](0,1) = 0.;
    ddphi[8](0,2) = -4.;
    ddphi[8](1,0) = 0.;
    ddphi[8](1,1) = 0.;
    ddphi[8](1,2) = -4.;
    ddphi[8](2,0) = -4.;
    ddphi[8](2,1) = -4.;
    ddphi[8](2,2) = -8.;

    ddphi[9](0,0) = 0.;
    ddphi[9](0,1) = -4.;
    ddphi[9](0,2) = 0.;
    ddphi[9](1,0) = -4.;
    ddphi[9](1,1) = -8.;
    ddphi[9](1,2) = -4.;
    ddphi[9](2,0) = 0.;
    ddphi[9](2,1) = -4.;
    ddphi[9](2,2) = 0.;

}

void ShapeTetrahedronQua::getCoordinates(MatrixDouble &coord) {

    coord(0,0) = 0.0;
    coord(1,0) = 0.0;
    coord(2,0) = 0.0;

    coord(0,1) = 0.0;
    coord(1,1) = 1.0;
    coord(2,1) = 0.0;

    coord(0,2) = 0.0;
    coord(1,2) = 0.0;
    coord(2,2) = 1.0;

    coord(0,3) = 1.0;
    coord(1,3) = 0.0;
    coord(2,3) = 0.0;

    coord(0,4) = 0.0;
    coord(1,4) = 0.5;
    coord(2,4) = 0.0;

    coord(0,5) = 0.0;
    coord(1,5) = 0.5;
    coord(2,5) = 0.5;

    coord(0,6) = 0.0;
    coord(1,6) = 0.0;
    coord(2,6) = 0.5;

    coord(0,7) = 0.5;
    coord(1,7) = 0.0;
    coord(2,7) = 0.0;

    coord(0,8) = 0.5;
    coord(1,8) = 0.5;
    coord(2,8) = 0.0;

    coord(0,9) = 0.5;
    coord(1,9) = 0.0;
    coord(2,9) = 0.5;

}
