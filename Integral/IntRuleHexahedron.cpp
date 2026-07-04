#include <stdio.h>
#include "IntRuleHexahedron.h"
#include "IntRule1d.h"

static int ComputingSymmetricCubatureRule(int order, MatrixDouble &Points, VecDouble &Weights);
void TransformBarycentricCoordInCartesianCoord(long double baryvec[], long double weightvec[], MatrixDouble &Points, VecDouble &Weight);

// Gauss quadrature points and weights for 1D integration
void gaussQuadrature(int numPoints, VecDouble& points, VecDouble& weights) {
    points.resize(numPoints);
    weights.resize(numPoints);

    // Weights and points for Gaussian quadrature with 3-point rule
    if (numPoints == 1) {
        points(0) = 0.0;
        weights(0) = 2.0;
    } else if (numPoints == 2) {
        points << -0.5773502691896257, 0.5773502691896257;
        weights << 1.0, 1.0;
    } else if (numPoints == 3) {
        points << -0.7745966692414834, 0.0, 0.7745966692414834;
        weights << 0.5555555555555556, 0.8888888888888888, 0.5555555555555556;
    } else if (numPoints == 4) {
        points << -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526;
        weights << 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538;
    } else if (numPoints == 5) {
        points << -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640;
        weights << 0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891;
    }
}

// Function to compute Gauss quadrature points and weights for hexahedron integration
std::pair<MatrixDouble, VecDouble> gaussQuadratureHexahedron(int numPointsPerAxis) {
    int numPoints = numPointsPerAxis * numPointsPerAxis * numPointsPerAxis;
    VecDouble gaussPoints1D, gaussWeights1D;
    gaussQuadrature(numPointsPerAxis, gaussPoints1D, gaussWeights1D);

    MatrixDouble points(3, numPoints);
    VecDouble weights(numPoints);

    int index = 0;
    for (int k = 0; k < numPointsPerAxis; ++k) {
        for (int j = 0; j < numPointsPerAxis; ++j) {
            for (int i = 0; i < numPointsPerAxis; ++i) {
                points.col(index) = Vector3R(gaussPoints1D(i), gaussPoints1D(j), gaussPoints1D(k));
                weights(index) = gaussWeights1D(i) * gaussWeights1D(j) * gaussWeights1D(k);
                ++index;
            }
        }
    }
    return std::make_pair(points, weights);
}



IntRuleHexahedron::IntRuleHexahedron() {
}

IntRuleHexahedron::IntRuleHexahedron(int order) : IntRule(order) {
    SetOrder(order);
}

void IntRuleHexahedron::SetOrder(int order) {

    fOrder = order;

#ifdef DEBUG_BUILD
    if (order < 0 || order > MaxOrder()) {
        PanicButton();
    }
#endif

    int nPoints = 2*order-1;
    if (order == 0) {
        nPoints = 1;
    }
    // fPoints.resize(nPoints*nPoints*nPoints,3);
    // fWeights.resize(nPoints);//is resized in GaulegQuad
    auto vals= gaussQuadratureHexahedron(nPoints);
    fPoints = vals.first.transpose();
    fWeights = vals.second;
    // gaussQuadrature3D(nPoints, fPoints, fWeights);

    // VecDouble coordAux(nPoints*nPoints*nPoints);
    // gaulegHexa(-1,1,coordAux,fWeights);

    // for (int i = 0; i < nPoints*nPoints*nPoints; i++){
    //     fPoints(i,0) = coordAux[i];
    //     fPoints(i,1) = coordAux[i+nPoints*nPoints*nPoints];
    //     fPoints(i,2) = coordAux[i+2*nPoints*nPoints*nPoints];
    // }
    // std::cout << "fPoints" << fPoints << std::endl;
    // std::cout << "fWeights" << fWeights << std::endl;
}

void IntRuleHexahedron::gaulegHexa(const REAL x1, const REAL x2, VecDouble &co, VecDouble &w) {
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

