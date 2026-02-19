#include "HierarquicalTriangle.h"
#include "HierarquicalOneD.h"

const int HierarquicalTriangle::Dimension;
const int HierarquicalTriangle::NCornerNodes;
const int HierarquicalTriangle::NSides;
const ElementType HierarquicalTriangle::ElType;

int HierarquicalTriangle::NSideNodes(int side) {
    if(side < 0 || side > 6) PanicButton();
    static int nsidenodes[7] = {1, 1, 1, 2, 2, 2, 3};
    return nsidenodes[side];
}

int HierarquicalTriangle::SideNodeLocIndex(int side, int node) {
    if(side<3 && node == 0) return side;
    if(side>=3 && side<6 && node <2) return (side-3+node) %3;
    if(side==6 && node <3) return node;
    std::cout << "HierarquicalTriangle::SideNodeIndex inconsistent side or node " << side
            << ' ' << node << std::endl;
    PanicButton();
    return -1;
}

/// returns the total number of shape functions
int HierarquicalTriangle::NShapeFunctions(VecInt &orders) {
    PanicButton();
    int nsf_tot = 0;
    for (int is=0; is<NSides; is++) {
        nsf_tot += NShapeFunctions(is, orders[is]);
    }
    
    return nsf_tot;
}

/// returns the number of shape functions associated with a side
int HierarquicalTriangle::NShapeFunctions(int side, int order){

  
    PanicButton();
    return -1;    
}


void HierarquicalTriangle::Shape(VecDouble &xi, VecDouble &phi, MatrixDouble &dphi, int order) {
    double xsi1 = xi[0];
    double xsi2 = xi[1];
    double xsi3 = 1. - xsi1 - xsi2;

    phi[0] = xsi3;
    phi[1] = xsi1;
    phi[2] = xsi2;   

    dphi(0,0) = -1.;
    dphi(1,0) = -1.;
    dphi(0,1) = 1.;
    dphi(1,1) = 0.;
    dphi(0,2) = 0.;
    dphi(1,2) = 1.;

}

void HierarquicalTriangle::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    PanicButton();
    return;
}

void HierarquicalTriangle::Chebyshev(double x, int num,MatrixDouble &phi,MatrixDouble &dphi){
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



void HierarquicalTriangle::Legendre(double x, int num,MatrixDouble &phi,MatrixDouble &dphi){
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
