#include "HierarquicalOneD.h"

const int HierarquicalOneD::Dimension;
const int HierarquicalOneD::NCornerNodes;
const int HierarquicalOneD::NSides;
const ElementType HierarquicalOneD::ElType;


int HierarquicalOneD::NSideNodes(int side) {
    if(side < 0 || side > 2) PanicButton();
    static int nsidenodes[3] = {1, 1, 2};
    return nsidenodes[side];
}

int HierarquicalOneD::SideNodeLocIndex(int side, int node) {
    if (side < 2 && node == 0) return side;
    if (side == 2 && node < 2) return node;
    std::cout << "HierarquicalOneD::SideNodeIndex inconsistent side or node " << side
            << ' ' << node << std::endl;
    PanicButton();
    return -1;
}



void HierarquicalOneD::Shape(VecDouble &xi, VecDouble &phi) {
    phi[0] = (1 - xi[0]) / 2.;
    phi[1] = (1 + xi[0]) / 2.;
}

void HierarquicalOneD::ShapeGradient(VecDouble &xi, MatrixDouble &dphi) {
    dphi(0,0) = -0.5;
    dphi(0,1) =  0.5;
}

void HierarquicalOneD::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    PanicButton();
    return;
}

void Chebyshev(double x, int num,MatrixDouble &phi,MatrixDouble &dphi){
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
