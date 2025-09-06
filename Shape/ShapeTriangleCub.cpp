#include "ShapeTriangleCub.h"

const int ShapeTriangleCub::Dimension;
const int ShapeTriangleCub::NSides;
const int ShapeTriangleCub::NShape;
const int ShapeTriangleCub::NCornerNodes;
const ElementType ShapeTriangleCub::ElType;

void ShapeTriangleCub::Shape(VecDouble &xi, VecDouble &phi) {

    double xsi1 = xi[0];
    double xsi2 = xi[1];
    double xsi3 = 1. - xsi1 - xsi2;

    phi[0] = (xsi1 * (3.0 * xsi1 - 2.0) * (3.0 * xsi1 - 1.0)) / 2.0;
    phi[1] = (xsi2 * (3.0 * xsi2 - 2.0) * (3.0 * xsi2 - 1.0)) / 2.0;
    phi[2] = -((xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0) * (3.0 * xsi2 + 3.0 * xsi1 - 1.0)) / 2.0;
    phi[3] = (9.0 * xsi1 * xsi2 * (3.0 * xsi1 - 1.0)) / 2.0;
    phi[4] = (9.0 * xsi1 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
    phi[5] = -(9.0 * xsi2 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 - 1.0)) / 2.0;
    phi[6] = (9.0 * xsi2 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0)) / 2.0;
    phi[7] = (9.0 * xsi1 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0)) / 2.0;
    phi[8] = -(9.0 * xsi1 * (3.0 * xsi1 - 1.0) * (xsi2 + xsi1 - 1.0)) / 2.0;
    phi[9] = -27.0 * xsi1 * xsi2 * (xsi2 + xsi1 - 1.0);
    // element conectivity
    //     2
    //     76
    //     895
    //     0341
        
}

void ShapeTriangleCub::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    const double xsi1 = xi[0];
    const double xsi2 = xi[1];
    const double xsi3 = 1. - xsi1 - xsi2;

    dphi(0,0) = (27.0 * xsi1 * xsi1 - 18.0 * xsi1 + 2.0) / 2.0;
    dphi(0,1) = 0.0;
    dphi(0,2) = -(27.0 * xsi2 * xsi2 + 54.0 * xsi1 * xsi2 - 36.0 * xsi2 + 27.0 * xsi1 * xsi1 - 36.0 * xsi1 + 11.0) / 2.0;
    dphi(0,3) = (9.0 * xsi2 * (6.0 * xsi1 - 1.0)) / 2.0;
    dphi(0,4) = (9.0 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
    dphi(0,5) = -(9.0 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
    dphi(0,6) = (9.0 * xsi2 * (6.0 * xsi2 + 6.0 * xsi1 - 5.0)) / 2.0;
    dphi(0,7) = (9.0 * (3.0 * xsi2 * xsi2 + 12.0 * xsi1 * xsi2 - 5.0 * xsi2 + 9.0 * xsi1 * xsi1 - 10.0 * xsi1 + 2.0)) / 2.0;
    dphi(0,8) = -(9.0 * (6.0 * xsi1 * xsi2 - xsi2 + 9.0 * xsi1 * xsi1 - 8.0 * xsi1 + 1.0)) / 2.0;
    dphi(0,9) = -27.0 * xsi2 * (xsi2 + 2.0 * xsi1 - 1.0);
    dphi(1,0) = 0.0;
    dphi(1,1) = (27.0 * xsi2 * xsi2 - 18.0 * xsi2 + 2) / 2.0;
    dphi(1,2) = -(27.0 * xsi2 * xsi2 + 54.0 * xsi1 * xsi2 - 36.0 * xsi2 + 27.0 * xsi1 * xsi1 - 36.0 * xsi1 + 11.0) / 2.0;
    dphi(1,3) = (9.0 * xsi1 * (3.0 * xsi1 - 1.0)) / 2.0;
    dphi(1,4) = (9.0 * xsi1 * (6.0 * xsi2 - 1.0)) / 2.0;
    dphi(1,5) = -(9.0 * (9.0 * xsi2 * xsi2 + 6.0 * xsi1 * xsi2 - 8.0 * xsi2 - xsi1 + 1.0)) / 2.0;
    dphi(1,6) = (9.0 * (9.0 * xsi2 * xsi2 + 12.0 * xsi1 * xsi2 - 10.0 * xsi2 + 3.0 * xsi1 * xsi1 - 5.0 * xsi1 + 2.0)) / 2.0;
    dphi(1,7) = (9.0 * xsi1 * (6.0 * xsi2 + 6.0 * xsi1 - 5.0)) / 2.0;
    dphi(1,8) = -(9.0 * xsi1 * (3.0 * xsi1 - 1.0)) / 2.0;
    dphi(1,9) = -27.0 * xsi1 * (2.0 * xsi2 + xsi1 - 1.0);
    // element conectivity
    //     2
    //     76
    //     895
    //     0341
}

void ShapeTriangleCub::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    const double xsi1 = xi[0];
    const double xsi2 = xi[1];
    
    ddphi[0](0,0) = -9. + 27. * xsi1;
    ddphi[0](0,1) = 0.;
    ddphi[0](1,0) = 0.;
    ddphi[0](1,1) = 0.;

    ddphi[1](0,0) = 0.;
    ddphi[1](0,1) = 0.;
    ddphi[1](1,0) = 0.;
    ddphi[1](1,1) = -9. + 27. * xsi2;

    ddphi[2](0,0) = 18. - 27.*xsi1 - 27.*xsi2;
    ddphi[2](0,1) = 18. - 27.*xsi1 - 27.*xsi2;
    ddphi[2](1,0) = 18. - 27.*xsi1 - 27.*xsi2;
    ddphi[2](1,1) = 18. - 27.*xsi1 - 27.*xsi2;

    ddphi[3](0,0) = 27. * xsi2;
    ddphi[3](0,1) = -4.5 + 27.*xsi1;
    ddphi[3](1,0) = -4.5 + 27.*xsi1;
    ddphi[3](1,1) = 0.;

    ddphi[4](0,0) = 0.;
    ddphi[4](0,1) = -4.5 + 27.*xsi2;
    ddphi[4](1,0) = -4.5 + 27.*xsi2;
    ddphi[4](1,1) = 27.*xsi1;

    ddphi[5](0,0) = 0.;
    ddphi[5](0,1) = 4.5 - 27.*xsi2;
    ddphi[5](1,0) = 4.5 - 27.*xsi2;
    ddphi[5](1,1) = 36. - 27.*xsi1 - 81.*xsi2;

    ddphi[6](0,0) = 27. * xsi2;
    ddphi[6](0,1) = -22.5 + 27.*xsi1+ 54.*xsi2;
    ddphi[6](1,0) = -22.5 + 27.*xsi1+ 54.*xsi2;
    ddphi[6](1,1) = -45. + 54. * xsi1 + 81. * xsi2;

    ddphi[7](0,0) = -45. + 81. * xsi1 + 54. * xsi2;
    ddphi[7](0,1) = -22.5 + 54.*xsi1+ 27.*xsi2;
    ddphi[7](1,0) = -22.5 + 54.*xsi1+ 27.*xsi2;
    ddphi[7](1,1) = 27. * xsi1;

    ddphi[8](0,0) = 36. - 27.*xsi2 - 81.*xsi1;
    ddphi[8](0,1) = 4.5 - 27.*xsi1;
    ddphi[8](1,0) = 4.5 - 27.*xsi1;
    ddphi[8](1,1) = 0.;

    ddphi[9](0,0) = -54. * xsi2;
    ddphi[9](0,1) = 27. - 54. * xsi1 - 54. * xsi2;
    ddphi[9](1,0) = 27. - 54. * xsi1 - 54. * xsi2;
    ddphi[9](1,1) = 54. * xsi1;
    //     2
    //     76
    //     895
    //     0341
    
}


void ShapeTriangleCub::getCoordinates(MatrixDouble &coord) {
    //     1
    //     54
    //     693
    //     2780
    coord(0,2) = 0.0;
    coord(1,2) = 0.0;

    coord(0,0) = 1.0;
    coord(1,0) = 0.0;

    coord(0,1) = 0.0;
    coord(1,1) = 1.0;

    coord(0,7) = 1./3.;
    coord(1,7) = 0.0;

    coord(0,8) = 2./3.;
    coord(1,8) = 0.0;

    coord(0,3) = 2./3.;
    coord(1,3) = 1./3.;

    coord(0,4) = 1./3.;
    coord(1,4) = 2./3.;

    coord(0,5) = 0.0;
    coord(1,5) = 2./3.;

    coord(0,6) = 0.0;
    coord(1,6) = 1./3.;

    coord(0,9) = 1./3.;
    coord(1,9) = 1./3.;
}
