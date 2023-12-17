#include <stdio.h>
#include "IntRuleHexahedron.h"
#include "IntRule1d.h"

static int ComputingSymmetricCubatureRule(int order, MatrixDouble &Points, VecDouble &Weights);
void TransformBarycentricCoordInCartesianCoord(long double baryvec[], long double weightvec[], MatrixDouble &Points, VecDouble &Weight);

IntRuleHexahedron::IntRuleHexahedron() {
}

IntRuleHexahedron::IntRuleHexahedron(int order) : IntRule(order) {
    SetOrder(order);
}

void IntRuleHexahedron::SetOrder(int order) {

    order = 2;
    fOrder = order;

    if (order < 0 || order > MaxOrder()) {
        PanicButton();
    }

    int nPoints = 2*order-1;
    if (order == 0) {
        nPoints = 1;
    }
    fPoints.resize(nPoints*nPoints*nPoints,3);
    fWeights.resize(nPoints);//is resized in GaulegQuad
    
    VecDouble coordAux(nPoints*nPoints*nPoints);
    gaulegHexa(-1,1,coordAux,fWeights);

    for (int i = 0; i < nPoints*nPoints*nPoints; i++){
        fPoints(i,0) = coordAux[i];
        fPoints(i,1) = coordAux[i+nPoints*nPoints*nPoints];
        fPoints(i,2) = coordAux[i+2*nPoints*nPoints*nPoints];
    }
    // std::cout << "fPoints" << fPoints << std::endl;
    // std::cout << "fWeights" << fWeights << std::endl;
}

void IntRuleHexahedron::gaulegHexa(const double x1, const double x2, VecDouble &co, VecDouble &w) {
    IntRule1d x;
    IntRule1d y;
    IntRule1d z;
    
    int n = w.size();   

    VecDouble cox(n);
    VecDouble coy(n);
    VecDouble coz(n);
    VecDouble wx(n);
    VecDouble wy(n);
    VecDouble wz(n);

    x.gauleg(x1, x2, cox, wx);
    y.gauleg(x1, x2, coy, wy);
    z.gauleg(x1, x2, coz, wz);
    
    co.resize(3*n*n*n);
    w.resize(n * n * n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                co[k + j * n + i * n * n] = cox[i];
                co[k + j * n + i * n * n + n*n*n] = coy[j];
                co[k + j * n + i * n * n + 2*n*n*n] = coz[k];
                w[n * n * i + n * j + k] = wx[i] * wy[j] * wz[k];
            }
        }
    }
}