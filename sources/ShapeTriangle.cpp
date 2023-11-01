#include "ShapeTriangle.h"

const int ShapeTriangle::Dimension;

void ShapeTriangle::Shape(VecDouble &xi, VecDouble &phi) const {

    double xsi1 = xi[0];
    double xsi2 = xi[1];
    double xsi3 = 1. - xsi1 - xsi2;
    switch (fOrder)
    {
    case 1:
        phi[0] = xsi3;
        phi[1] = xsi1;
        phi[2] = xsi2;
        // element conectivity
        //     2
        //     01
        break;
    case 2:
        phi[0] = xsi3 * (2.0 * xsi3 - 1.0);
        phi[1] = xsi1 * (2.0 * xsi1 - 1.0);
        phi[2] = xsi2 * (2.0 * xsi2 - 1.0);
        phi[3] = 4.0 * xsi3 * xsi1;
        phi[4] = 4.0 * xsi1 * xsi2;
        phi[5] = 4.0 * xsi2 * xsi3;
        // element conectivity
        //     2
        //     54
        //     031
        break;

    case 3:
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
        break;
    
    default:
        PanicButton();
        break;
    }
}

void ShapeTriangle::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) const {
    const double xsi1 = xi[0];
    const double xsi2 = xi[1];
    const double xsi3 = 1. - xsi1 - xsi2;

    switch (fOrder)
    {
    case 1:
        dphi(0,0) = -1.;
        dphi(0,1) = -1.;
        dphi(1,0) = 1.;
        dphi(1,1) = 0.;
        dphi(2,0) = 0.;
        dphi(2,1) = 1.;
        // element conectivity
        //     2
        //     01
        break;
    case 2:
        dphi(1,0) = 4. * xsi1 - 1.;
        dphi(1,1) = 0.;
        dphi(2,0) = 0.;
        dphi(2,1) = 4. * xsi2 - 1.;
        dphi(0,0) = -4. * xsi3 + 1.;
        dphi(0,1) = -4. * xsi3 + 1.;
        dphi(4,0) = 4. * xsi2;
        dphi(4,1) = 4. * xsi1;
        dphi(5,0) = -4. * xsi2;
        dphi(5,1) = 4. * (xsi3 - xsi2);
        dphi(3,0) = 4. * (xsi3 - xsi1);
        dphi(3,1) = -4. * xsi1;
        // element conectivity
        //     2
        //     54
        //     031
        break;
    case 3:
        dphi(0,0) = (27.0 * xsi1 * xsi1 - 18.0 * xsi1 + 2.0) / 2.0;
        dphi(1,0) = 0.0;
        dphi(2,0) = -(27.0 * xsi2 * xsi2 + 54.0 * xsi1 * xsi2 - 36.0 * xsi2 + 27.0 * xsi1 * xsi1 - 36.0 * xsi1 + 11.0) / 2.0;
        dphi(3,0) = (9.0 * xsi2 * (6.0 * xsi1 - 1.0)) / 2.0;
        dphi(4,0) = (9.0 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
        dphi(5,0) = -(9.0 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
        dphi(6,0) = (9.0 * xsi2 * (6.0 * xsi2 + 6.0 * xsi1 - 5.0)) / 2.0;
        dphi(7,0) = (9.0 * (3.0 * xsi2 * xsi2 + 12.0 * xsi1 * xsi2 - 5.0 * xsi2 + 9.0 * xsi1 * xsi1 - 10.0 * xsi1 + 2.0)) / 2.0;
        dphi(8,0) = -(9.0 * (6.0 * xsi1 * xsi2 - xsi2 + 9.0 * xsi1 * xsi1 - 8.0 * xsi1 + 1.0)) / 2.0;
        dphi(9,0) = -27.0 * xsi2 * (xsi2 + 2.0 * xsi1 - 1.0);
        dphi(0,1) = 0.0;
        dphi(1,1) = (27.0 * xsi2 * xsi2 - 18.0 * xsi2 + 2) / 2.0;
        dphi(2,1) = -(27.0 * xsi2 * xsi2 + 54.0 * xsi1 * xsi2 - 36.0 * xsi2 + 27.0 * xsi1 * xsi1 - 36.0 * xsi1 + 11.0) / 2.0;
        dphi(3,1) = (9.0 * xsi1 * (3.0 * xsi1 - 1.0)) / 2.0;
        dphi(4,1) = (9.0 * xsi1 * (6.0 * xsi2 - 1.0)) / 2.0;
        dphi(5,1) = -(9.0 * (9.0 * xsi2 * xsi2 + 6.0 * xsi1 * xsi2 - 8.0 * xsi2 - xsi1 + 1.0)) / 2.0;
        dphi(6,1) = (9.0 * (9.0 * xsi2 * xsi2 + 12.0 * xsi1 * xsi2 - 10.0 * xsi2 + 3.0 * xsi1 * xsi1 - 5.0 * xsi1 + 2.0)) / 2.0;
        dphi(7,1) = (9.0 * xsi1 * (6.0 * xsi2 + 6.0 * xsi1 - 5.0)) / 2.0;
        dphi(8,1) = -(9.0 * xsi1 * (3.0 * xsi1 - 1.0)) / 2.0;
        dphi(9,1) = -27.0 * xsi1 * (2.0 * xsi2 + xsi1 - 1.0);
        // element conectivity
        //     2
        //     76
        //     895
        //     0341
        break;

    default:
        PanicButton();
        break;
    }
}

void ShapeTriangle::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) const {
    const double xsi1 = xi[0];
    const double xsi2 = xi[1];
    switch (fOrder)
    {
    case 1:
        return;
        break;
    case 2:
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
        break;
    case 3:
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
        break;
    
    default:
        PanicButton();
        break;
    }

}


void ShapeTriangle::getCoordinates(MatrixDouble &coord) const {

    switch (fOrder)
    {
    case 1:
        coord(0,0) = 0.0;
        coord(1,0) = 0.0;

        coord(0,1) = 1.0;
        coord(1,1) = 0.0;

        coord(0,2) = 0.0;
        coord(1,2) = 1.0;
        break;
    case 2:
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
        break;
    case 3:
        std::cout << "Not implemented\n";
        PanicButton();
        break;
    
    default:
        PanicButton();
        break;
    }

}
