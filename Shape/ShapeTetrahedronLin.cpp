#include "ShapeTetrahedronLin.h"

const int ShapeTetrahedronLin::Dimension;
const int ShapeTetrahedronLin::NSides;
const int ShapeTetrahedronLin::NShape;
const int ShapeTetrahedronLin::NCornerNodes;
const int ShapeTetrahedronLin::NEdges;
const int ShapeTetrahedronLin::NFaces;
const int ShapeTetrahedronLin::NVolumes;
const ElementType ShapeTetrahedronLin::ElType;

void ShapeTetrahedronLin::Shape(VecDouble &xi, VecDouble &phi, MatrixDouble &dphi, VecInt orders) {

    REAL xsi1 = xi[0];
    REAL xsi2 = xi[1];
    REAL xsi3 = xi[2];
    REAL xsi4 = 1.0 - xsi1 - xsi2 - xsi3;

    phi[0] = 1.0 - xsi1 - xsi2 - xsi3;
    phi[1] = xsi1;
    phi[2] = xsi2;
    phi[3] = xsi3;
    
    dphi(0,0) = -1.0;
    dphi(1,0) = -1.0;
    dphi(2,0) = -1.0;
    dphi(0,1) = 1.0;
    dphi(1,1) = 0.0;
    dphi(2,1) = 0.0;
    dphi(0,2) = 0.0;
    dphi(1,2) = 1.0;
    dphi(2,2) = 0.0;
    dphi(0,3) = 0.0;
    dphi(1,3) = 0.0;
    dphi(2,3) = 1.0;
}

void ShapeTetrahedronLin::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    return;

}

void ShapeTetrahedronLin::getCoordinates(MatrixDouble &coord) {

    coord(0,0) = 0.0;
    coord(1,0) = 0.0;
    coord(2,0) = 0.0;

    coord(0,1) = 1.0;
    coord(1,1) = 0.0;
    coord(2,1) = 0.0;

    coord(0,2) = 0.0;
    coord(1,2) = 1.0;
    coord(2,2) = 0.0;

    coord(0,3) = 0.0;
    coord(1,3) = 0.0;
    coord(2,3) = 1.0;
}
