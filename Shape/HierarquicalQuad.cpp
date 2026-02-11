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


void HierarquicalQuad::Shape(VecDouble &xi, VecDouble &phi, int order) {
    
    if(order <= 0) PanicButton();

    MatrixDouble phi0(order+1,1),phi1(order+1,1);
    MatrixDouble dphi0(1,order+1),dphi1(1,order+1);

    HierarquicalOneD::Legendre(xi[0],order,phi0,dphi0);
    HierarquicalOneD::Legendre(xi[1],order,phi1,dphi1);

    for (int i=0;i<order;i++) {
        for (int j=0;j<order;j++) {
            int index = i*order+j;
            phi(index,0) =  phi0(i,0) * phi1(j,0);
            // dphi(0,index) = dphi0(0,i)* phi1(j,0);
            // dphi(1,index) =  phi0(i,0)*dphi1(0,j);
        }
    }

}

void HierarquicalQuad::ShapeGradient(VecDouble &xi, MatrixDouble &dphi, int order) {
    dphi(0,0) = -0.5;
    dphi(0,1) =  0.5;

    if (order >= 2) {
        MatrixDouble phih(order+1,1);
        MatrixDouble dphih(1,order+1);
        Legendre(xi[0],order,phih,dphih);
        for(int i=2;i<=order;i++) dphi(0,i) = dphih(0,i)-dphih(0,i-2);
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
