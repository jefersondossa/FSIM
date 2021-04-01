//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2020 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------------BOUNDARY QUADRATURE POINTS---------------------------
//------------------------------------------------------------------------------

#ifndef BOUND_INTEG_QUADRATURE_H
#define BOUND_INTEG_QUADRATURE_H

#include "BoundaryShapeFunction.hpp"

/// Defines the quadrature rule for the boundary integration
template<int DIM, int DEG>
class BoundaryIntegQuadrature{
private:
    double pi = M_PI;  //Pi

public:

    /// Returns the index of the first integration point
    /// @return first integration point index
    double* begin() {
        return std::begin(pointWeight);
    }

    /// Returns the index of the last integration point
    /// @return last integration point index
    double* end() {
        return std::end(pointWeight);
    }

    /// Returns the integration point coordinate
    /// @param int integration point index @param int adimensional direction
    /// @return integration point adimensional coordinates
    double PointList(int i, int j) {return pointCoord[i][j];}; 

    /// Retuns the integration point weight
    /// @param int integration point index @return integration point weight
    double WeightList(int i) {return pointWeight[i];};

    /// Sets the integration quadrature
    void setQuadrature();

    /// Boundary integration quadrature constructor
    BoundaryIntegQuadrature(){
        setQuadrature();
    }


private:
    ///List of integration points coordinates
    double pointCoord[-3*DIM-8*DEG+4*DIM*DEG+8][DIM-1];

    ///List of integration points weights
    double pointWeight[-3*DIM-8*DEG+4*DIM*DEG+8];
};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------
template<>
void BoundaryIntegQuadrature<2,1>::setQuadrature(){

    // 2 Gauss points
    int numPoints = 2;
    double xmga, xlga, zga, p1ga, p2ga, p3ga, ppga, z1ga;
    int mga;
    int nga = numPoints;

    mga=(nga+1.)/2.;
    xmga=0.0;
    xlga=1.0;
    
    for (int iga=1; iga<=mga; iga++) {
        zga = std::cos(pi*(double(iga)-0.25)/(double(nga)+0.5));
    g1:
        p1ga = 1.0;
        p2ga = 0.0;
        for (int jga=1; jga <= nga; jga++) {
            p3ga = p2ga;
            p2ga = p1ga;
            p1ga = ((2.0*double(jga)-1.0)*zga*p2ga-(double(jga)-1.0)*p3ga)/(double(jga));
        };
     
        ppga = nga*(zga*p1ga-p2ga)/(zga*zga-1.0);
        z1ga = zga;
        zga = z1ga-p1ga/ppga;
        
        if (std::fabs(zga-z1ga) > 1.0e-15) goto g1;
        
        pointCoord[iga-1][0] = xmga-xlga*zga;
        pointCoord[nga-iga][0] = xmga + xlga*zga;
        pointWeight[iga-1] = 2.0*xlga/((1.0-zga*zga)*ppga*ppga);
        pointWeight[nga-iga] = pointWeight[iga-1];
    };
    return;
}; 

template<>
void BoundaryIntegQuadrature<2,2>::setQuadrature(){

    // 2 Gauss points
    int numPoints = 2;
    double xmga, xlga, zga, p1ga, p2ga, p3ga, ppga, z1ga;
    int mga;
    int nga = numPoints;

    mga=(nga+1.)/2.;
    xmga=0.0;
    xlga=1.0;
    
    for (int iga=1; iga<=mga; iga++) {
        zga = std::cos(pi*(double(iga)-0.25)/(double(nga)+0.5));
    g1:
        p1ga = 1.0;
        p2ga = 0.0;
        for (int jga=1; jga <= nga; jga++) {
            p3ga = p2ga;
            p2ga = p1ga;
            p1ga = ((2.0*double(jga)-1.0)*zga*p2ga-(double(jga)-1.0)*p3ga)/(double(jga));
        };
     
        ppga = nga*(zga*p1ga-p2ga)/(zga*zga-1.0);
        z1ga = zga;
        zga = z1ga-p1ga/ppga;
        
        if (std::fabs(zga-z1ga) > 1.0e-15) goto g1;
        
        pointCoord[iga-1][0] = xmga-xlga*zga;
        pointCoord[nga-iga][0] = xmga + xlga*zga;
        pointWeight[iga-1] = 2.0*xlga/((1.0-zga*zga)*ppga*ppga);
        pointWeight[nga-iga] = pointWeight[iga-1];
    };
    return;
}; 


template<>
void BoundaryIntegQuadrature<3,1>::setQuadrature(){

    // 3 Hammer Points
    pointCoord[0][0] = 0.166666666666666666666666666666667e+0;
    pointCoord[0][1] = 0.666666666666666666666666666666667e+0;

    pointCoord[1][0] = 0.666666666666666666666666666666667e+0;
    pointCoord[1][1] = 0.166666666666666666666666666666667e+0;

    pointCoord[2][0] = 0.166666666666666666666666666666667e+0;
    pointCoord[2][1] = 0.166666666666666666666666666666667e+0;

    pointWeight[0]  =  0.166666666666666666666666666666667e+0;
    pointWeight[1]  =  0.166666666666666666666666666666667e+0;
    pointWeight[2]  =  0.166666666666666666666666666666667e+0;

    return;
}

template<>
void BoundaryIntegQuadrature<3,2>::setQuadrature(){

    // 7 Hammer Points
    double g1 = (6. - std::sqrt(15.)) / 21.;
    double g2 = (6. + std::sqrt(15.)) / 21.;

    double w1 = (155. - std::sqrt(15.)) / 2400.;
    double w2 = (155. + std::sqrt(15.)) / 2400.;

    pointCoord[0][0] = 1. / 3.;
    pointCoord[0][1] = 1. / 3.;
        
    pointCoord[1][0] = 1. - 2.*g1;
    pointCoord[1][1] = g1;
      
    pointCoord[2][0] = g1;
    pointCoord[2][1] = 1. - 2.*g1;
      
    pointCoord[3][0] = g1;
    pointCoord[3][1] = g1;
      
    pointCoord[4][0] = g2;
    pointCoord[4][1] = g2;
      
    pointCoord[5][0] = 1. - 2.*g2;
    pointCoord[5][1] = g2;
      
    pointCoord[6][0] = g2;
    pointCoord[6][1] = 1. - 2.*g2;

    pointWeight[0] = 0.11250;
    pointWeight[1] = w1;
    pointWeight[2] = w1;
    pointWeight[3] = w1;
    pointWeight[4] = w2;
    pointWeight[5] = w2;
    pointWeight[6] = w2;

    return;
}

#endif
