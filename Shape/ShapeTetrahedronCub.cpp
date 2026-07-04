#include "ShapeTetrahedronCub.h"

const int ShapeTetrahedronCub::Dimension;
const int ShapeTetrahedronCub::NSides;
const int ShapeTetrahedronCub::NShape;
const int ShapeTetrahedronCub::NCornerNodes;
const int ShapeTetrahedronCub::NEdges;
const int ShapeTetrahedronCub::NFaces;
const int ShapeTetrahedronCub::NVolumes;
const ElementType ShapeTetrahedronCub::ElType;

void ShapeTetrahedronCub::Shape(VecDouble &xi, VecDouble &phi, MatrixDouble &dphi, VecInt orders) {

    REAL xsi1 = xi[0];
    REAL xsi2 = xi[1];
    REAL xsi3 = xi[2];
    REAL xsi4 = 1.0 - xsi1 - xsi2 - xsi3;

    phi[0] = (1.0 / 2.0) * (3.0 * xsi4 - 1.0) * (3.0 * xsi4 - 2.0) * xsi4;
    phi[1] = (1.0 / 2.0) * (3.0 * xsi1 - 1.0) * (3.0 * xsi1 - 2.0) * xsi1;
    phi[2] = (1.0 / 2.0) * (3.0 * xsi2 - 1.0) * (3.0 * xsi2 - 2.0) * xsi2;
    phi[3] = (1.0 / 2.0) * (3.0 * xsi3 - 1.0) * (3.0 * xsi3 - 2.0) * xsi3;
    phi[4] = (9.0 / 2.0) * (3.0 * xsi4 - 1.0) * xsi1 * xsi4;
    phi[5] = (9.0 / 2.0) * (3.0 * xsi1 - 1.0) * xsi1 * xsi4;
    phi[6] = (9.0 / 2.0) * (3.0 * xsi1 - 1.0) * xsi1 * xsi2;
    phi[7] = (9.0 / 2.0) * (3.0 * xsi2 - 1.0) * xsi1 * xsi2;
    phi[8] = (9.0 / 2.0) * (3.0 * xsi2 - 1.0) * xsi2 * xsi4;
    phi[9] = (9.0 / 2.0) * (3.0 * xsi4 - 1.0) * xsi2 * xsi4;
    phi[10] = (9.0 / 2.0) * (3.0 * xsi3 - 1.0) * xsi3 * xsi4;
    phi[11] = (9.0 / 2.0) * (3.0 * xsi4 - 1.0) * xsi3 * xsi4;
    phi[12] = (9.0 / 2.0) * (3.0 * xsi3 - 1.0) * xsi2 * xsi3;
    phi[13] = (9.0 / 2.0) * (3.0 * xsi2 - 1.0) * xsi2 * xsi3;
    phi[14] = (9.0 / 2.0) * (3.0 * xsi3 - 1.0) * xsi1 * xsi3;
    phi[15] = (9.0 / 2.0) * (3.0 * xsi1 - 1.0) * xsi1 * xsi3;
    phi[16] = 27.0 * xsi1 * xsi2 * xsi4;
    phi[17] = 27.0 * xsi1 * xsi3 * xsi4;
    phi[18] = 27.0 * xsi2 * xsi3 * xsi4;
    phi[19] = 27.0 * xsi1 * xsi2 * xsi3;

    dphi(0,0 ) = 0.50 * (-11.0 + 36.0 * xsi1 - 27.0 * xsi1 * xsi1 + 36.0 * xsi2 - 54.0 * xsi1 * xsi2 - 27.0 * xsi2 * xsi2 + 36.0 * xsi3 - 54.0 * xsi1 * xsi3 - 54.0 * xsi2 * xsi3 - 27.0 * xsi3 * xsi3);
    dphi(0,1 ) = 0.50 * (2.0 - 18.0 * xsi1 + 27.0 * xsi1 * xsi1);
    dphi(0,2 ) = 0.0;
    dphi(0,3 ) = 0.0;
    dphi(0,4 ) = 9.0 / 2.0 * (2.0 - 10.0 * xsi1 + 9.0 * xsi1 * xsi1 - 5.0 * xsi2 + 12.0 * xsi1 * xsi2 + 3.0 * xsi2 * xsi2 - 5.0 * xsi3 + 12.0 * xsi1 * xsi3 + 6.0 * xsi2 * xsi3 + 3.0 * xsi3 * xsi3);
    dphi(0,5 ) = -9.0 / 2.0 * (1.0 - 8.0 * xsi1 + 9.0 * xsi1 * xsi1 - xsi2 + 6.0 * xsi1 * xsi2 - xsi3 + 6.0 * xsi1 * xsi3);
    dphi(0,6 ) = 9.0 / 2.0 * (-1.0 + 6.0 * xsi1) * xsi2;
    dphi(0,7 ) = 9.0 / 2.0 * xsi2 * (-1.0 + 3.0 * xsi2);
    dphi(0,8 ) = -9.0 / 2.0 * xsi2 * (-1.0 + 3.0 * xsi2);
    dphi(0,9 ) = 9.0 / 2.0 * xsi2 * (-5.0 + 6.0 * xsi1 + 6.0 * xsi2 + 6.0 * xsi3);
    dphi(0,10) = -9.0 / 2.0 * xsi3 * (-1.0 + 3.0 * xsi3);
    dphi(0,11) = 9.0 / 2.0 * xsi3 * (-5.0 + 6.0 * xsi1 + 6.0 * xsi2 + 6.0 * xsi3);
    dphi(0,12) = 0.0;
    dphi(0,13) = 0.0;
    dphi(0,14) = 9.0 / 2.0 * xsi3 * (-1.0 + 3.0 * xsi3);
    dphi(0,15) = 9.0 / 2.0 * (-1.0 + 6.0 * xsi1) * xsi3;
    dphi(0,16) = -27.0 * xsi2 * (-1.0 + 2.0 * xsi1 + xsi2 + xsi3);
    dphi(0,17) = -27.0 * xsi3 * (-1.0 + 2.0 * xsi1 + xsi2 + xsi3);
    dphi(0,18) = -27.0 * xsi2 * xsi3;
    dphi(0,19) = 27.0 * xsi2 * xsi3;

    dphi(1,0 ) = 0.50 * (-11.0 + 36.0 * xsi1 - 27.0 * xsi1 * xsi1 + 36.0 * xsi2 - 54.0 * xsi1 * xsi2 - 27.0 * xsi2 * xsi2 + 36.0 * xsi3 - 54.0 * xsi1 * xsi3 - 54.0 * xsi2 * xsi3 - 27.0 * xsi3 * xsi3);
    dphi(1,1 ) = 0.0;
    dphi(1,2 ) = 0.50 * (2.0 - 18.0 * xsi2 + 27.0 * xsi2 * xsi2);
    dphi(1,3 ) = 0.0;
    dphi(1,4 ) = 9.0 / 2.0 * xsi1 * (-5.0 + 6.0 * xsi1 + 6.0 * xsi2 + 6.0 * xsi3);
    dphi(1,5 ) = -9.0 / 2.0 * xsi1 * (-1.0 + 3.0 * xsi1);
    dphi(1,6 ) = 9.0 / 2.0 * xsi1 * (-1.0 + 3.0 * xsi1);
    dphi(1,7 ) = 9.0 / 2.0 * xsi1 * (-1.0 + 6.0 * xsi2);
    dphi(1,8 ) = -9.0 / 2.0 * (1.0 - xsi1 - 8.0 * xsi2 - xsi3 + 6.0 * xsi1 * xsi2 + 6.0 * xsi2 * xsi3 + 9.0 * xsi2 * xsi2);
    dphi(1,9 ) = 9.0 / 2.0 * (2.0 - 5.0 * xsi1 + 3.0 * xsi1 * xsi1 - 10.0 * xsi2 + 12.0 * xsi1 * xsi2 + 9.0 * xsi2 * xsi2 - 5.0 * xsi3 + 6.0 * xsi1 * xsi3 + 12.0 * xsi2 * xsi3 + 3.0 * xsi3 * xsi3);
    dphi(1,10) = -9.0 / 2.0 * xsi3 * (-1.0 + 3.0 * xsi3);
    dphi(1,11) = 9.0 / 2.0 * xsi3 * (-5.0 + 6.0 * xsi1 + 6.0 * xsi2 + 6.0 * xsi3);
    dphi(1,12) = 9.0 / 2.0 * (-1.0 + 3.0 * xsi3) * xsi3;
    dphi(1,13) = 9.0 / 2.0 * (-1.0 + 6.0 * xsi2) * xsi3;
    dphi(1,14) = 0.0;
    dphi(1,15) = 0.0;
    dphi(1,16) = -27.0 * xsi1 * (-1.0 + xsi1 + 2.0 * xsi2 + xsi3);
    dphi(1,17) = -27.0 * xsi1 * xsi3;
    dphi(1,18) = -27.0 * xsi3 * (-1.0 + xsi1 + 2.0 * xsi2 + xsi3);
    dphi(1,19) = 27.0 * xsi1 * xsi3;

    dphi(2,0 ) = 0.50 * (-11.0 + 36.0 * xsi1 - 27.0 * xsi1 * xsi1 + 36.0 * xsi2 - 54.0 * xsi1 * xsi2 - 27.0 * xsi2 * xsi2 + 36.0 * xsi3 - 54.0 * xsi1 * xsi3 - 54.0 * xsi2 * xsi3 - 27.0 * xsi3 * xsi3);
    dphi(2,1 ) = 0.0;
    dphi(2,2 ) = 0.0;
    dphi(2,3 ) = 0.50 * (2.0 - 18.0 * xsi3 + 27.0 * xsi3 * xsi3);
    dphi(2,4 ) = 9.0 / 2.0 * xsi1 * (-5.0 + 6.0 * xsi1 + 6.0 * xsi2 + 6.0 * xsi3);
    dphi(2,5 ) = -9.0 / 2.0 * xsi1 * (-1.0 + 3.0 * xsi1);
    dphi(2,6 ) = 0.0;
    dphi(2,7 ) = 0.0;
    dphi(2,8 ) = -9.0 / 2.0 * xsi2 * (-1.0 + 3.0 * xsi2);
    dphi(2,9 ) = 9.0 / 2.0 * xsi2 * (-5.0 + 6.0 * xsi1 + 6.0 * xsi2 + 6.0 * xsi3);
    dphi(2,10) = -9.0 / 2.0 * (1.0 - xsi1 - xsi2 - 8.0 * xsi3 + 6.0 * xsi1 * xsi3 + 6.0 * xsi2 * xsi3 + 9.0 * xsi3 * xsi3);
    dphi(2,11) = 9.0 / 2.0 * (2.0 - 5.0 * xsi1 + 3.0 * xsi1 * xsi1 - 5.0 * xsi2 + 6.0 * xsi1 * xsi2 + 3.0 * xsi2 * xsi2 - 10.0 * xsi3 + 12.0 * xsi1 * xsi3 + 12.0 * xsi2 * xsi3 + 9.0 * xsi3 * xsi3);
    dphi(2,12) = 9.0 / 2.0 * xsi2 * (-1.0 + 6.0 * xsi3);
    dphi(2,13) = 9.0 / 2.0 * xsi2 * (-1.0 + 3.0 * xsi2);
    dphi(2,14) = 9.0 / 2.0 * xsi1 * (-1.0 + 6.0 * xsi3);
    dphi(2,15) = 9.0 / 2.0 * xsi1 * (-1.0 + 3.0 * xsi1);
    dphi(2,16) = -27.0 * xsi1 * xsi2;
    dphi(2,17) = -27.0 * xsi1 * (-1.0 + xsi1 + 2.0 * xsi3 + xsi2);
    dphi(2,18) = -27.0 * xsi2 * (-1.0 + xsi1 + 2.0 * xsi3 + xsi2);
    dphi(2,19) = 27.0 * xsi1 * xsi2;
    
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

void ShapeTetrahedronCub::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    const REAL xsi1 = xi[0];
    const REAL xsi2 = xi[1];
    const REAL xsi3 = xi[2];

    ddphi[0](0,0) = 18. - 27.*xsi1 - 27.*xsi2 - 27.*xsi3;
    ddphi[0](0,1) = 18. - 27.*xsi1 - 27.*xsi2 - 27.*xsi3;
    ddphi[0](0,2) = 18. - 27.*xsi1 - 27.*xsi2 - 27.*xsi3;
    ddphi[0](1,0) = 18. - 27.*xsi1 - 27.*xsi2 - 27.*xsi3;
    ddphi[0](1,1) = 18. - 27.*xsi1 - 27.*xsi2 - 27.*xsi3;
    ddphi[0](1,2) = 18. - 27.*xsi1 - 27.*xsi2 - 27.*xsi3;
    ddphi[0](2,0) = 18. - 27.*xsi1 - 27.*xsi2 - 27.*xsi3;
    ddphi[0](2,1) = 18. - 27.*xsi1 - 27.*xsi2 - 27.*xsi3;
    ddphi[0](2,2) = 18. - 27.*xsi1 - 27.*xsi2 - 27.*xsi3;
    
    ddphi[1](0,0) = -9. + 27. * xsi1;
    ddphi[1](0,1) = 0.;
    ddphi[1](0,2) = 0.;
    ddphi[1](1,0) = 0.;
    ddphi[1](1,1) = 0.;
    ddphi[1](1,2) = 0.;
    ddphi[1](2,0) = 0.;
    ddphi[1](2,1) = 0.;
    ddphi[1](2,2) = 0.;

    ddphi[2](0,0) = 0.;
    ddphi[2](0,1) = 0.;
    ddphi[2](0,2) = 0.;
    ddphi[2](1,0) = 0.;
    ddphi[2](1,1) = -9. + 27. * xsi2;
    ddphi[2](1,2) = 0.;
    ddphi[2](2,0) = 0.;
    ddphi[2](2,1) = 0.;
    ddphi[2](2,2) = 0.;

    ddphi[3](0,0) = 0.;
    ddphi[3](0,1) = 0.;
    ddphi[3](0,2) = 0.;
    ddphi[3](1,0) = 0.;
    ddphi[3](1,1) = 0.;
    ddphi[3](1,2) = 0.;
    ddphi[3](2,0) = 0.;
    ddphi[3](2,1) = 0.;
    ddphi[3](2,2) = -9. + 27. * xsi3;

    ddphi[4](0,0) = -45. + 81.* xsi1 + 54. * xsi2 + 54. * xsi3;
    ddphi[4](0,1) = -22.5 + 54. * xsi1 + 27. * xsi2 + 27. * xsi3;
    ddphi[4](0,2) = -22.5 + 54. * xsi1 + 27. * xsi2 + 27. * xsi3;
    ddphi[4](1,0) = -22.5 + 54. * xsi1 + 27. * xsi2 + 27. * xsi3;
    ddphi[4](1,1) = 27. * xsi1;
    ddphi[4](1,2) = 27. * xsi1;
    ddphi[4](2,0) = -22.5 + 54. * xsi1 + 27. * xsi2 + 27. * xsi3;
    ddphi[4](2,1) = 27. * xsi1;
    ddphi[4](2,2) = 27. * xsi1;

    ddphi[5](0,0) = 36. - 81. * xsi1 - 27. * xsi2 - 27. * xsi3;
    ddphi[5](0,1) = 4.5 - 27. * xsi1;
    ddphi[5](0,2) = 4.5 - 27. * xsi1;
    ddphi[5](1,0) = 4.5 - 27. * xsi1;
    ddphi[5](1,1) = 0.;
    ddphi[5](1,2) = 0.;
    ddphi[5](2,0) = 4.5 - 27. * xsi1;
    ddphi[5](2,1) = 0.;
    ddphi[5](2,2) = 0.;

    ddphi[6](0,0) = 27. * xsi2;
    ddphi[6](0,1) = -4.5 + 27. * xsi1;
    ddphi[6](0,2) = 0.;
    ddphi[6](1,0) = -4.5 + 27. * xsi1;
    ddphi[6](1,1) = 0.;
    ddphi[6](1,2) = 0.;
    ddphi[6](2,0) = 0.;
    ddphi[6](2,1) = 0.;
    ddphi[6](2,2) = 0.;

    ddphi[7](0,0) = 0.;
    ddphi[7](0,1) = -4.5 + 27. * xsi2;
    ddphi[7](0,2) = 0.;
    ddphi[7](1,0) = -4.5 + 27. * xsi2;
    ddphi[7](1,1) = 27. * xsi1;
    ddphi[7](1,2) = 0.;
    ddphi[7](2,0) = 0.;
    ddphi[7](2,1) = 0.;
    ddphi[7](2,2) = 0.;

    ddphi[8](0,0) = 0.;
    ddphi[8](0,1) = 4.5 - 27. * xsi2;
    ddphi[8](0,2) = 0.;
    ddphi[8](1,0) = 4.5 - 27. * xsi2;
    ddphi[8](1,1) = -27. * xsi1 - 81. * xsi2 -27.* xsi3;
    ddphi[8](1,2) = 4.5 - 27. * xsi2;
    ddphi[8](2,0) = 0.;
    ddphi[8](2,1) = 4.5 - 27. * xsi2;
    ddphi[8](2,2) = 0.;

    ddphi[9](0,0) = 27. * xsi2;
    ddphi[9](0,1) = -22.5 + 27. * xsi1 + 54. * xsi2 + 27. * xsi3;
    ddphi[9](0,2) = 27. * xsi2;
    ddphi[9](1,0) = -22.5 + 27. * xsi1 + 54. * xsi2 + 27. * xsi3;
    ddphi[9](1,1) = -45. + 54. * xsi1 + 81. * xsi2 + 54. * xsi3;
    ddphi[9](1,2) = -22.5 + 27. * xsi1 + 54. * xsi2 + 27. * xsi3;
    ddphi[9](2,0) = 27. * xsi2;
    ddphi[9](2,1) = -22.5 + 27. * xsi1 + 54. * xsi2 + 27. * xsi3;
    ddphi[9](2,2) = 27. * xsi2;

    ddphi[10](0,0) = 0.;
    ddphi[10](0,1) = 0.;
    ddphi[10](0,2) = 4.5 - 27. * xsi3;
    ddphi[10](1,0) = 0.;
    ddphi[10](1,1) = 0.;
    ddphi[10](1,2) = 4.5 - 27. * xsi3;
    ddphi[10](2,0) = 4.5 - 27. * xsi3;
    ddphi[10](2,1) = 4.5 - 27. * xsi3;
    ddphi[10](2,2) = 36. - 27. * xsi1 - 27. * xsi2 - 81. * xsi3;

    ddphi[11](0,0) = 27. * xsi3;
    ddphi[11](0,1) = 27. * xsi3;
    ddphi[11](0,2) = -22.5 + 27. * xsi1 + 27. * xsi2 + 54. * xsi3;
    ddphi[11](1,0) = 27. * xsi3;
    ddphi[11](1,1) = 27. * xsi3;
    ddphi[11](1,2) = -22.5 + 27. * xsi1 + 27. * xsi2 + 54. * xsi3;
    ddphi[11](2,0) = -22.5 + 27. * xsi1 + 27. * xsi2 + 54. * xsi3;
    ddphi[11](2,1) = -22.5 + 27. * xsi1 + 27. * xsi2 + 54. * xsi3;
    ddphi[11](2,2) = -45. + 54. * xsi1 + 54. * xsi2 + 81. * xsi3;

    ddphi[12](0,0) = 0.;
    ddphi[12](0,1) = 0.;
    ddphi[12](0,2) = 0.;
    ddphi[12](1,0) = 0.;
    ddphi[12](1,1) = 0.;
    ddphi[12](1,2) = -4.5 + 27. * xsi3;
    ddphi[12](2,0) = 0.;
    ddphi[12](2,1) = -4.5 + 27. * xsi3;
    ddphi[12](2,2) = 27. * xsi2;

    ddphi[13](0,0) = 0.;
    ddphi[13](0,1) = 0.;
    ddphi[13](0,2) = 0.;
    ddphi[13](1,0) = 0.;
    ddphi[13](1,1) = 27. * xsi3;
    ddphi[13](1,2) = -4.5 + 27. * xsi2;
    ddphi[13](2,0) = 0.;
    ddphi[13](2,1) = -4.5 + 27. * xsi2;
    ddphi[13](2,2) = 0.;

    ddphi[14](0,0) = 0.;
    ddphi[14](0,1) = 0.;
    ddphi[14](0,2) = -4.5 + 27. * xsi3;
    ddphi[14](1,0) = 0.;
    ddphi[14](1,1) = 0.;
    ddphi[14](1,2) = 0.;
    ddphi[14](2,0) = -4.5 + 27. * xsi3;
    ddphi[14](2,1) = 0.;
    ddphi[14](2,2) = 27. * xsi1;

    ddphi[15](0,0) = 27. * xsi3;
    ddphi[15](0,1) = 0.;
    ddphi[15](0,2) = -4.5 + 27. * xsi1;
    ddphi[15](1,0) = 0.;
    ddphi[15](1,1) = 0.;
    ddphi[15](1,2) = 0.;
    ddphi[15](2,0) = -4.5 + 27. * xsi1;
    ddphi[15](2,1) = 0.;
    ddphi[15](2,2) = 0.;

    ddphi[16](0,0) = -54. * xsi2;
    ddphi[16](0,1) = 27. - 54. * xsi1 - 54. * xsi2 - 27. * xsi3;
    ddphi[16](0,2) = -27. * xsi2;
    ddphi[16](1,0) = 27. - 54. * xsi1 - 54. * xsi2 - 27. * xsi3;
    ddphi[16](1,1) = -54. * xsi1;
    ddphi[16](1,2) = -27. * xsi1;
    ddphi[16](2,0) = -27. * xsi2;
    ddphi[16](2,1) = -27. * xsi1;
    ddphi[16](2,2) = 0.;

    ddphi[17](0,0) = -54. * xsi3;
    ddphi[17](0,1) = -27. * xsi3;
    ddphi[17](0,2) = 27. - 54. * xsi1 - 27. * xsi2 - 54. * xsi3;
    ddphi[17](1,0) = -27. * xsi3;
    ddphi[17](1,1) = 0.;
    ddphi[17](1,2) = -27. * xsi1;
    ddphi[17](2,0) = 27. - 54. * xsi1 - 27. * xsi2 - 54. * xsi3;
    ddphi[17](2,1) = -27. * xsi1;
    ddphi[17](2,2) = -54. * xsi1;

    ddphi[18](0,0) = 0.;
    ddphi[18](0,1) = -27. * xsi3;
    ddphi[18](0,2) = -27. * xsi2;
    ddphi[18](1,0) = -27. * xsi3;
    ddphi[18](1,1) = -54. * xsi3;
    ddphi[18](1,2) = 27. - 27. * xsi1 - 54. * xsi2 - 54. * xsi3;
    ddphi[18](2,0) = -27. * xsi2;
    ddphi[18](2,1) = 27. - 27. * xsi1 - 54. * xsi2 - 54. * xsi3;
    ddphi[18](2,2) = -54. * xsi2;

    ddphi[19](0,0) = 0.;
    ddphi[19](0,1) = 27. * xsi3;
    ddphi[19](0,2) = 27. * xsi2;
    ddphi[19](1,0) = 27. * xsi3;
    ddphi[19](1,1) = 0.;
    ddphi[19](1,2) = 27. * xsi1;
    ddphi[19](2,0) = 27. * xsi2;
    ddphi[19](2,1) = 27. * xsi1;
    ddphi[19](2,2) = 0.;
           
}

void ShapeTetrahedronCub::getCoordinates(MatrixDouble &coord) {

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

    coord(0,4) = 1./3.;
    coord(1,4) = 0.0;
    coord(2,4) = 0.0;

    coord(0,5) = 2./3.;
    coord(1,5) = 0.0;
    coord(2,5) = 0.0;

    coord(0,6) = 0.0;
    coord(1,6) = 1./3.;
    coord(2,6) = 0.0;

    coord(0,7) = 0.0;
    coord(1,7) = 2./3.;
    coord(2,7) = 0.0;

    coord(0,8) = 0.0;
    coord(1,8) = 0.0;
    coord(2,8) = 1./3.;

    coord(0,9) = 0.0;
    coord(1,9) = 0.0;
    coord(2,9) = 2./3.;

    coord(0,10) = 2./3.;
    coord(1,10) = 1./3.;
    coord(2,10) = 0.0;

    coord(0,11) = 1./3.;
    coord(1,11) = 2./3.;
    coord(2,11) = 0.0;

    coord(0,12) = 0.0;
    coord(1,12) = 2./3.;
    coord(2,12) = 1./3.;

    coord(0,13) = 0.0;
    coord(1,13) = 1./3.;
    coord(2,13) = 2./3.;

    coord(0,14) = 2./3.;
    coord(1,14) = 0.0;
    coord(2,14) = 1./3.;

    coord(0,15) = 1./3.;
    coord(1,15) = 0.0;
    coord(2,15) = 2./3.;

    coord(0,16) = 1./3.;
    coord(1,16) = 1./3.;
    coord(2,16) = 0.0;

    coord(0,17) = 1./3.;
    coord(1,17) = 0.0;
    coord(2,17) = 1./3.;

    coord(0,18) = 0.0;
    coord(1,18) = 1./3.;
    coord(2,18) = 1./3.;

    coord(0,19) = 1./3.;
    coord(1,19) = 1./3.;
    coord(2,19) = 1./3.;

}
