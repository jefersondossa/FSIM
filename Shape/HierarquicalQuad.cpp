#include "HierarquicalQuad.h"
#include "HierarquicalOneD.h"

const int HierarquicalQuad::Dimension;
const int HierarquicalQuad::NCornerNodes;
const int HierarquicalQuad::NSides;
const ElementType HierarquicalQuad::ElType;

int HierarquicalQuad::NSideNodes(int side) {
    if(side < 0 || side > 8) PanicButton();
    static int nsidenodes[9] = {1, 1, 1, 1, 2, 2, 2, 2, 4};
    return nsidenodes[side];
}

int HierarquicalQuad::SideNodeLocIndex(int side, int node) {
    if (side < 4 && node == 0) return side;
    if(side >=4 && side < 8 && node < 2) return (side+node)%4;
    if(side == 8 && node < 4) return node;
    std::cout << "HierarquicalQuad::SideNodeIndex inconsistent side or node " << side
            << ' ' << node << std::endl;
    PanicButton();
    return -1;
}

/// returns the total number of shape functions
int HierarquicalQuad::NShapeFunctions(VecInt &orders) {
    
    int nsf_tot = 0;
    for (int is=0; is<NSides; is++) {
        nsf_tot += NShapeFunctions(is, orders[is]);
    }
    
    return nsf_tot;
}

/// returns the number of shape functions associated with a side
int HierarquicalQuad::NShapeFunctions(int side, int order){

    if(order < 1 ) PanicButton();
    
    if(side<4) return 1;//0 a 4

    if(side<8) return (order-1);//4 a 7
    if(side==8) {
        return ((order-1)*(order-1));
    }


    std::cout << "Shape1d::NShapeFunctions : Wrong side " << side << "\n";
    PanicButton();
    return -1;    
}


void HierarquicalQuad::Shape(VecDouble &xi, VecDouble &phi, MatrixDouble &dphi, int order) {
    
    if(order <= 0) PanicButton();
    if(order >=3){
        std::cout << "HierarquicalQuad::Shape : order not implemented, the side (edge) shape functions need to be compatibilized " << order << "\n";
        PanicButton();
    }

    phi[0] = 0.25 * (1 - xi[0])*(1 - xi[1]);
    phi[1] = 0.25 * (1 + xi[0])*(1 - xi[1]);
    phi[2] = 0.25 * (1 + xi[0])*(1 + xi[1]);
    phi[3] = 0.25 * (1 - xi[0])*(1 + xi[1]);

    dphi(0,0) = -0.25 * (1. - xi[1]);
    dphi(1,0) = -0.25 * (1. - xi[0]);
    dphi(0,1) =  0.25 * (1. - xi[1]);
    dphi(1,1) = -0.25 * (1. + xi[0]);
    dphi(0,2) =  0.25 * (1. + xi[1]);
    dphi(1,2) =  0.25 * (1. + xi[0]);
    dphi(0,3) = -0.25 * (1. + xi[1]);
    dphi(1,3) =  0.25 * (1. - xi[0]);

    if (order == 1) return;

    int nshape1d = HierarquicalOneD::NShapeFunctions(order);
    VecDouble phi0(nshape1d), phi1(nshape1d);
    MatrixDouble dphi0(1,nshape1d), dphi1(1,nshape1d);
    VecDouble xi0(1), xi1(1);
    xi0[0] = xi[0];
    xi1[0] = xi[1];
    HierarquicalOneD::Shape(xi0, phi0, dphi0, order);
    HierarquicalOneD::Shape(xi1, phi1, dphi1, order);

    int count = 4;
    // Side 4
    for(int j=2; j < nshape1d; j++) {
        phi[count] = phi0[j]*phi1[0];
        dphi(0,count) = dphi0(0,j)*phi1[0];
        dphi(1,count) = phi0[j]*dphi1(0,0);
        count++;
    }
    // Side 5
    for(int j=2; j < nshape1d; j++) {
        phi[count] = phi0[1]*phi1[j];
        dphi(0,count) = dphi0(0,1)*phi1[j];
        dphi(1,count) = phi0[1]*dphi1(0,j);
        count++;
    }
    // Side 6
    for(int j=2; j < nshape1d; j++) {
        phi[count] = phi0[j]*phi1[1];
        dphi(0,count) = dphi0(0,j)*phi1[1];
        dphi(1,count) = phi0[j]*dphi1(0,1);
        count++;
    }
    // Side 7
    for(int j=2; j < nshape1d; j++) {
        phi[count] = phi0[0]*phi1[j];
        dphi(0,count) = dphi0(0,0)*phi1[j];
        dphi(1,count) = phi0[0]*dphi1(0,j);
        count++;
    }
    // Side 8
    for(int j=2; j < nshape1d; j++) {
        for(int k=2; k < nshape1d; k++) {
            phi[count] = phi0[j]*phi1[k];
            dphi(0,count) = dphi0(0,j)*phi1[k];
            dphi(1,count) = phi0[j]*dphi1(0,k);
            count++;
        }
    }    
}

void HierarquicalQuad::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    PanicButton();
    return;
}

void HierarquicalQuad::Chebyshev(double x, int num,MatrixDouble &phi,MatrixDouble &dphi){
    // Quadratic or higher shape functions
    if(num <= 0) return;
    phi.setZero();
    dphi.setZero();
    phi(0,0) = 1.0;
    dphi(0,0) = 0.0;
    if(num == 1) return;
    phi(1,0) = x;
    dphi(0,1) = 1.0;
    int ord;
    // dphi.Print("DphisAntes = ",std::cout,EMathematicaInput);

    for(ord = 2;ord<num;ord++) {
        phi(ord,0) = 2.0*x*phi(ord-1,0) - phi(ord-2,0);
        dphi(0,ord) = 2.0*x*dphi(0,ord-1) + 2.0*phi(ord-1,0) - dphi(0,ord-2);
    }
    // dphi.Print("DphisDepois = ",std::cout,EMathematicaInput);
}



void HierarquicalQuad::Legendre(double x, int num,MatrixDouble &phi,MatrixDouble &dphi){
    // Quadratic or higher shape functions
    if(num <= 0) return;
    phi.setZero();
    dphi.setZero();
    phi(0,0) = 1.0;
    dphi(0,0) = 0.0;
    if(num == 1) return;
    phi(1,0) = x;
    dphi(0,1) = 1.0;
    int ord;
    // dphi.Print("DphisAntes = ",std::cout,EMathematicaInput);

    for(ord = 2;ord<=num;ord++) {
        phi(ord,0) = ((2.0*ord-1.0)*x*phi(ord-1,0) - (ord-1.0)*phi(ord-2,0))/ord;
        dphi(0,ord) = dphi(0,ord-2) + (2.0*ord-1.0)*phi(ord-1,0);
    }
    // dphi.Print("DphisDepois = ",std::cout,EMathematicaInput);
}
