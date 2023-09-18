#include "BoundaryIntegrationQuadrature.h"


//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------
void BoundaryIntegQuadrature::setQuadrature(){

    if (pointCoord.cols() == 2){
        int numPoints = 0;
        switch (fOrder)
        {
        case 1:
        case 2:
            numPoints = 2;
            break;
        case 3:
            numPoints = 3;
            break;
        default:
            break;
        }
        //Gauss points
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
            
            pointCoord(iga-1,0) = xmga-xlga*zga;
            pointCoord(nga-iga,0) = xmga + xlga*zga;
            pointWeight[iga-1] = 2.0*xlga/((1.0-zga*zga)*ppga*ppga);
            pointWeight[nga-iga] = pointWeight[iga-1];
        };

    } else if (pointCoord.cols() == 3){
        switch (fOrder)
        {
        case 1:
            pointCoord(0,0) = 0.166666666666666666666666666666667e+0;
            pointCoord(0,1) = 0.666666666666666666666666666666667e+0;

            pointCoord(1,0) = 0.666666666666666666666666666666667e+0;
            pointCoord(1,1) = 0.166666666666666666666666666666667e+0;

            pointCoord(2,0) = 0.166666666666666666666666666666667e+0;
            pointCoord(2,1) = 0.166666666666666666666666666666667e+0;

            pointWeight[0]  =  0.166666666666666666666666666666667e+0;
            pointWeight[1]  =  0.166666666666666666666666666666667e+0;
            pointWeight[2]  =  0.166666666666666666666666666666667e+0;
            break;
        case 2:
            {
                // 7 Hammer Points
                double g1 = (6. - std::sqrt(15.)) / 21.;
                double g2 = (6. + std::sqrt(15.)) / 21.;

                double w1 = (155. - std::sqrt(15.)) / 2400.;
                double w2 = (155. + std::sqrt(15.)) / 2400.;

                pointCoord(0,0) = 1. / 3.;
                pointCoord(0,1) = 1. / 3.;
                    
                pointCoord(1,0) = 1. - 2.*g1;
                pointCoord(1,1) = g1;
                
                pointCoord(2,0) = g1;
                pointCoord(2,1) = 1. - 2.*g1;
                
                pointCoord(3,0) = g1;
                pointCoord(3,1) = g1;
                
                pointCoord(4,0) = g2;
                pointCoord(4,1) = g2;
                
                pointCoord(5,0) = 1. - 2.*g2;
                pointCoord(5,1) = g2;
                
                pointCoord(6,0) = g2;
                pointCoord(6,1) = 1. - 2.*g2;

                pointWeight[0] = 0.11250;
                pointWeight[1] = w1;
                pointWeight[2] = w1;
                pointWeight[3] = w1;
                pointWeight[4] = w2;
                pointWeight[5] = w2;
                pointWeight[6] = w2;
            }
            break;
        case 3:
            std::cout << "Not Implemented\n";
            PanicButton();
            break;
        
        default:
            PanicButton();
            break;
        }

    } else {
        PanicButton();
    }


    return;
}; 

