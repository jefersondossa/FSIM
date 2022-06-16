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

#include "IntegrationQuadrature.h"
#include "IntegrationQuadrature11.h"
#include "Node.h"

/// Computes and stores the shape functions and its derivatives
template<int DIM, int DEG>
class DomainIntegration{
public:
    DomainIntegration(){
        IntegQuadrature<DIM,DEG>        quad;
        QuadShapeFunction<DIM,DEG>      shapeFunction;
        IntegQuadratureSpecial<DIM,DEG> Squad;
        int nElNodes = 3*(DIM*DEG-DEG)-2*DIM+4;

        int nIntegPointsNormal = quad.getNumberOfIntegrationPoints();
        int nIntegPointsSpecial = Squad.getNumberOfIntegrationPoints();

        phi_.resize(nElNodes,nIntegPointsNormal);
        dphi_.resize(nElNodes);
        // dphi_ = new double**[nElNodes];
        for (int i = 0; i < nElNodes; i++){ 
            dphi_[i].resize(DIM,nIntegPointsNormal);
            // dphi_[i] = new double*[DIM];
            // for (int j = 0; j < DIM; j++) dphi_[i][j] = new double[nIntegPointsNormal];
        }

        phiS_.resize(nElNodes,nIntegPointsSpecial);

        dphiS_ = new double**[nElNodes];
        for (int i = 0; i < nElNodes; i++){ 
            dphiS_[i] = new double*[DIM];
            for (int j = 0; j < DIM; j++) dphiS_[i][j] = new double[nIntegPointsSpecial];
        }

        VecDouble xsi(DIM);
        int index = 0;
        
        VecDouble phiAux_(nElNodes);

        MatrixDouble dphiAux(nElNodes,DIM);
        
        for(double* it = quad.begin(); it != quad.end(); it++){
            //Defines the integration points adimentional coordinates
            for (int i = 0; i < DIM; i++) xsi[i] = quad.PointList(index,i);       
            //Shape functions
            shapeFunction.evaluate(xsi,phiAux_);
            for (int i = 0; i < nElNodes; i++) phi_(i,index) = phiAux_[i];
            //Derivatives
            shapeFunction.evaluateGradient(xsi,dphiAux);
            for (int i = 0; i < nElNodes; i++)
                for (int j = 0; j < DIM; j++)
                    dphi_[i](j,index) = dphiAux(i,j);

            index++;
        }
        index = 0;
        for(double* it = Squad.begin(); it != Squad.end(); it++){
            //Defines the integration points adimentional coordinates
            for (int i = 0; i < DIM; i++) xsi[i] = Squad.PointList(index,i);       
            //Shape functions
            shapeFunction.evaluate(xsi,phiAux_);
            for (int i = 0; i < nElNodes; i++) phiS_(i,index) = phiAux_[i];
            //Derivatives
            shapeFunction.evaluateGradient(xsi,dphiAux);
            for (int i = 0; i < nElNodes; i++)
                for (int j = 0; j < DIM; j++)
                    dphiS_[i][j][index] = dphiAux(i,j);

            index++;
        }

    }

    MatrixDouble phi_;
    std::vector<MatrixDouble> dphi_;
    // double ***dphi_;

    MatrixDouble phiS_;
    double ***dphiS_;

private:

    
};






#endif
