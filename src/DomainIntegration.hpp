//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2020 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------DOMAIN INTEGRATION SHAPE FUNCTIONS AND DERIVATIVES--------------
//------------------------------------------------------------------------------

#ifndef DOMAIN_INTEGRATION_H
#define DOMAIN_INTEGRATION_H

#include "IntegrationQuadrature.hpp"
#include "Node.hpp"

/// Computes and stores the shape functions and its derivatives
template<int DIM, int DEG>
class DomainIntegration{
public:
    DomainIntegration(){
        IntegQuadrature<DIM,DEG>      quad;
        QuadShapeFunction<DIM,DEG>    shapeFunction;

        double xsi[DIM] = {};
        int index = 0;
            
        int nElNodes = 3*(DIM*DEG-DEG)-2*DIM+4;

        double phiAux_[nElNodes] = {};

        double **dphiAux;
        dphiAux = new double*[nElNodes];
        for (int i = 0; i < nElNodes; ++i) dphiAux[i] = new double[DIM];

        for(double* it = quad.begin(); it != quad.end(); it++){
            //Defines the integration points adimentional coordinates
            for (int i = 0; i < DIM; i++) xsi[i] = quad.PointList(index,i);       
            //Shape functions
            shapeFunction.evaluate(xsi,phiAux_);
            for (int i = 0; i < nElNodes; i++) phi_[i][index] = phiAux_[i];
            //Derivatives
            shapeFunction.evaluateGradient(xsi,dphiAux);
            for (int i = 0; i < nElNodes; i++)
                for (int j = 0; j < DIM; j++)
                    dphi_[i][j][index] = dphiAux[i][j];

            index++;
        }

        for (int i = 0; i < nElNodes; ++i) delete [] dphiAux[i];
        delete [] dphiAux;
    }

    double phi_[3*(DIM*DEG-DEG)-2*DIM+4][-5*DIM-8*DEG+6*DIM*DEG+9];
    double dphi_[3*(DIM*DEG-DEG)-2*DIM+4][DIM][-5*DIM-8*DEG+6*DIM*DEG+9];

private:

    
};






#endif
