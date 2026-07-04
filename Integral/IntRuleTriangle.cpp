#include <iostream> 
#include "IntRuleTriangle.h"

IntRuleTriangle::IntRuleTriangle(){

}

IntRuleTriangle::IntRuleTriangle(int order) {
    SetOrder(order);
}

void IntRuleTriangle::SetOrder(int order) {
    fOrder = order;

#ifdef DEBUG_BUILD
    if (order < 0 || order > MaxOrder()) {
        PanicButton();
    }
#endif

    switch (order)
    {
    case 0:
    case 1:
        fPoints.resize(1,2);
        fWeights.resize(1);
        fPoints(0,0) = 1./3.;
        fPoints(0,1) = 1./3.;
        fWeights[0] = 0.5;        
        break;
    case 2:
        fPoints.resize(4,2);
        fWeights.resize(4);
        fPoints(0,0) = 1./3.;
        fPoints(0,1) = 1./3.;
        fPoints(1,0) = 0.2;
        fPoints(1,1) = 0.6;
        fPoints(2,0) = 0.2;
        fPoints(2,1) = 0.2;
        fPoints(3,0) = 0.6;
        fPoints(3,1) = 0.2;
        fWeights[0]  =  -0.5625/2.;
        fWeights[1]  =  0.520833333333333333333333333333333/2.;
        fWeights[2]  =  0.520833333333333333333333333333333/2.;
        fWeights[3]  =  0.520833333333333333333333333333333/2.;
        break;
    case 3:
    case 4:
    case 5:
    case 6:
    {
        fPoints.resize(7,2);
        fWeights.resize(7);
        REAL g1 = (6. - std::sqrt(15.)) / 21.;
        REAL g2 = (6. + std::sqrt(15.)) / 21.;

        REAL w1 = (155. - std::sqrt(15.)) / 2400.;
        REAL w2 = (155. + std::sqrt(15.)) / 2400.;

        fPoints(0,0) = 1. / 3.;
        fPoints(0,1) = 1. / 3.;
            
        fPoints(1,0) = 1. - 2.*g1;
        fPoints(1,1) = g1;
        
        fPoints(2,0) = g1;
        fPoints(2,1) = 1. - 2.*g1;
        
        fPoints(3,0) = g1;
        fPoints(3,1) = g1;
        
        fPoints(4,0) = g2;
        fPoints(4,1) = g2;
        
        fPoints(5,0) = 1. - 2.*g2;
        fPoints(5,1) = g2;
        
        fPoints(6,0) = g2;
        fPoints(6,1) = 1. - 2.*g2;

        fWeights[0] = 0.11250;
        fWeights[1] = w1;
        fWeights[2] = w1;
        fWeights[3] = w1;
        fWeights[4] = w2;
        fWeights[5] = w2;
        fWeights[6] = w2;
    }
        break;
    case 7:
    case 8:
    case 9:
        {
        fPoints.resize(12,2);
        fWeights.resize(12);
        fPoints(0,0)     =  0.630890144915022283403316028708192e-1;
        fPoints(0,1)     =  0.873821971016995543319336794258362e+0;
        fPoints(1,0)     =  0.873821971016995543319336794258362e+0;
        fPoints(1,1)     =  0.630890144915022283403316028708192e-1;
        fPoints(2,0)     =  0.630890144915022283403316028708192e-1;
        fPoints(2,1)     =  0.630890144915022283403316028708192e-1;
        fPoints(3,0)     =  0.249286745170910421291638553107019e+0;
        fPoints(3,1)     =  0.501426509658179157416722893785962e+0;
        fPoints(4,0)     =  0.501426509658179157416722893785962e+0;
        fPoints(4,1)     =  0.249286745170910421291638553107019e+0;
        fPoints(5,0)     =  0.249286745170910421291638553107019e+0;
        fPoints(5,1)     =  0.249286745170910421291638553107019e+0;
        fPoints(6,0)     =  0.531450498448169473532496716313982e-1;
        fPoints(6,1)     =  0.636502499121398647230142594412050e+0;
        fPoints(7,0)     =  0.636502499121398647230142594412050e+0;
        fPoints(7,1)     =  0.531450498448169473532496716313982e-1;
        fPoints(8,0)     =  0.310352451033784405416607733956552e+0;
        fPoints(8,1)     =  0.636502499121398647230142594412050e+0;
        fPoints(9,0)     =  0.636502499121398647230142594412050e+0;
        fPoints(9,1)     =  0.310352451033784405416607733956552e+0;
        fPoints(10,0)    =  0.310352451033784405416607733956552e+0;
        fPoints(10,1)    =  0.531450498448169473532496716313982e-1;
        fPoints(11,0)    =  0.531450498448169473532496716313982e-1;
        fPoints(11,1)    =  0.310352451033784405416607733956552e+0;

        fWeights[0]   =  0.254224531851034084604684045534345e-1;
        fWeights[1]   =  0.254224531851034084604684045534345e-1;
        fWeights[2]   =  0.254224531851034084604684045534345e-1;
        fWeights[3]   =  0.583931378631896830126448056927897e-1;
        fWeights[4]   =  0.583931378631896830126448056927897e-1;
        fWeights[5]   =  0.583931378631896830126448056927897e-1;
        fWeights[6]   =  0.414255378091867875967767282102212e-1;
        fWeights[7]   =  0.414255378091867875967767282102212e-1;
        fWeights[8]   =  0.414255378091867875967767282102212e-1;
        fWeights[9]   =  0.414255378091867875967767282102212e-1;
        fWeights[10]  =  0.414255378091867875967767282102212e-1;
        fWeights[11]  =  0.414255378091867875967767282102212e-1;
        }
        break;
    default:
        PanicButton();
        break;
    }
}
