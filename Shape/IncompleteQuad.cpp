#include "IncompleteQuad.h"
#include "HierarquicalOneD.h"

const int IncompleteQuad::Dimension;
const int IncompleteQuad::NCornerNodes;
const int IncompleteQuad::NSides;
const ElementType IncompleteQuad::ElType;

int IncompleteQuad::NSideNodes(int side) {
    if(side < 0 || side > 8) PanicButton();
    static int nsidenodes[9] = {1, 1, 1, 1, 2, 2, 2, 2, 4};
    return nsidenodes[side];
}

int IncompleteQuad::SideNodeLocIndex(int side, int node) {
    if (side < 4 && node == 0) return side;
    if(side >=4 && side < 8 && node < 2) return (side+node)%4;
    if(side == 8 && node < 4) return node;
    std::cout << "IncompleteQuad::SideNodeIndex inconsistent side or node " << side
            << ' ' << node << std::endl;
    PanicButton();
    return -1;
}

/// returns the total number of shape functions
int IncompleteQuad::NShapeFunctions(VecInt &orders) {
    
    int nsf_tot = 0;
    for (int is=0; is<NSides; is++) {
        nsf_tot += NShapeFunctions(is, orders[is]);
    }
    
    return nsf_tot;
}

/// returns the number of shape functions associated with a side
int IncompleteQuad::NShapeFunctions(int side, int order){

    if(order < 0) PanicButton();
    
    if(side<3) return 1;//0 a 2
    if (side >= 3 && side < 9) return 0;

    std::cout << "IncompleteQuad::NShapeFunctions : Wrong side " << side << "\n";
    PanicButton();
    return -1;    
}


void IncompleteQuad::Shape(VecDouble &xi, VecDouble &phi, MatrixDouble &dphi, VecInt orders) {
    
    int minorder = orders.minCoeff();
    int maxorder = orders.maxCoeff();
    if(minorder < 0) PanicButton();

    if (minorder == 0) {
        phi[0] = 1.0;
        dphi.setZero();
        return;
    }
    
    if(maxorder >=2){
        std::cout << "IncompleteQuad::Shape : order not implemented, the side (edge) shape functions need to be compatibilized "  << "\n";
        PanicButton();
    }

    phi[0] = 1.;
    phi[1] = xi[0];
    phi[2] = xi[1];

    dphi(0,0) = 0.;
    dphi(1,0) = 0.;
    dphi(0,1) = 1.;
    dphi(1,1) = 0.;
    dphi(0,2) = 0.;
    dphi(1,2) = 1.;
    if (maxorder == 1) return;

}

void IncompleteQuad::ShapeHessian(VecDouble &xi, std::vector<MatrixDouble > &ddphi) {
    PanicButton();
    return;
}
