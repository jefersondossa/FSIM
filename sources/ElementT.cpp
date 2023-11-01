#include "ElementT.h"
// #include "Boundary.h"
// #include "CompMesh.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------CREATES AN AUXILIARY FINITE ELEMENT OF DIMENSION---------------
//----------------------DIM-1 FOR THE BOUNDARY INTEGRATION----------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::getBoundaryNodes(int *nodesb_){
    if (fMesh->Dimension() == 2){
        switch (fMesh->GetDefaultOrder())
        {
        case 1:
            if(fSideInBoundary == 0){
                nodesb_[0] = fConnect[2]; 
                nodesb_[1] = fConnect[1]; 
            }else{
                if(fSideInBoundary == 1){
                    nodesb_[0] = fConnect[0]; 
                    nodesb_[1] = fConnect[2]; 
                }else{
                    nodesb_[0] = fConnect[1];
                    nodesb_[1] = fConnect[0];
                };        
            };
            break;
        case 2:
            //!!!!!!!ATENÇÃO ESSAS CONECTIVIDADES FORAM ALTERADAS EM RELAÇÃO AO PROGRAMA Fluid.
            if(fSideInBoundary == 0){
                nodesb_[0] = fConnect[1]; 
                nodesb_[1] = fConnect[2]; 
                nodesb_[2] = fConnect[4];         
            }else{
                if(fSideInBoundary == 1){
                    nodesb_[0] = fConnect[2]; 
                    nodesb_[1] = fConnect[0]; 
                    nodesb_[2] = fConnect[5]; 
                }else{
                    nodesb_[0] = fConnect[0];
                    nodesb_[1] = fConnect[1];
                    nodesb_[2] = fConnect[3];
                };        
            };
            break;
        case 3:
            if(fSideInBoundary == 0){
                nodesb_[0] = fConnect[1]; 
                nodesb_[1] = fConnect[2]; 
                nodesb_[2] = fConnect[5];         
                nodesb_[3] = fConnect[6];         
            }else{
                if(fSideInBoundary == 1){
                    nodesb_[0] = fConnect[2]; 
                    nodesb_[1] = fConnect[0]; 
                    nodesb_[2] = fConnect[7]; 
                    nodesb_[3] = fConnect[8]; 
                }else{
                    nodesb_[0] = fConnect[0];
                    nodesb_[1] = fConnect[1];
                    nodesb_[2] = fConnect[3];
                    nodesb_[3] = fConnect[4];
                };        
            };
            break;
        default:
            PanicButton();
            break;
        }
    } else if (fMesh->Dimension() == 3){
        switch (fMesh->GetDefaultOrder())
        {
        case 1:
            if(fSideInBoundary == 0){
                nodesb_[0] = fConnect[1]; 
                nodesb_[1] = fConnect[2]; 
                nodesb_[2] = fConnect[3];         
            }else{
                if(fSideInBoundary == 1){
                    nodesb_[0] = fConnect[0]; 
                    nodesb_[1] = fConnect[3]; 
                    nodesb_[2] = fConnect[2]; 
                }else{
                    if(fSideInBoundary == 2){
                        nodesb_[0] = fConnect[0];
                        nodesb_[1] = fConnect[1];
                        nodesb_[2] = fConnect[3];
                    }else{
                        nodesb_[0] = fConnect[0];
                        nodesb_[1] = fConnect[2];
                        nodesb_[2] = fConnect[1];
                    }
                };        
            };
            break;
        case 2:
            if(fSideInBoundary == 0){
                nodesb_[0] = fConnect[2]; 
                nodesb_[1] = fConnect[3]; 
                nodesb_[2] = fConnect[1];
                nodesb_[3] = fConnect[9]; 
                nodesb_[4] = fConnect[8]; 
                nodesb_[5] = fConnect[5];
            }else{
                if(fSideInBoundary == 1){
                    nodesb_[0] = fConnect[0]; 
                    nodesb_[1] = fConnect[3]; 
                    nodesb_[2] = fConnect[2];
                    nodesb_[3] = fConnect[7]; 
                    nodesb_[4] = fConnect[9]; 
                    nodesb_[5] = fConnect[6];
                }else{
                    if(fSideInBoundary == 2){
                        nodesb_[0] = fConnect[3]; 
                        nodesb_[1] = fConnect[0]; 
                        nodesb_[2] = fConnect[1];
                        nodesb_[3] = fConnect[7]; 
                        nodesb_[4] = fConnect[4]; 
                        nodesb_[5] = fConnect[8];
                    }else{
                        nodesb_[0] = fConnect[0]; 
                        nodesb_[1] = fConnect[2]; 
                        nodesb_[2] = fConnect[1];
                        nodesb_[3] = fConnect[6]; 
                        nodesb_[4] = fConnect[5]; 
                        nodesb_[5] = fConnect[4];
                    }
                };        
            };
            break;
        case 3:
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
}
//------------------------------------------------------------------------------
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------


template<class tshape>
void ElementT<tshape>::ComputeIntPointDistFunction(VecDouble &nodalval) {
    
    int DIM = tshape::Dimension;
    VecDouble xsi(DIM);
    // ShapeFunction shapeQuad(DIM,DEG);
    // VecDouble phi_(fMesh->NElNodes());
    
    IntegQuadrature nQuad(DIM,DEG);
    // for(int i = 0; i < nQuad.getNumberOfIntegrationPoints(); i++) {
    //     intPointWeightFunctionPrev[i] = intPointWeightFunction[i];
    //     intPointWeightFunction[i] = 0.;
    // }

    int index=0;

    for(int it = 0; it < nQuad.getNumberOfIntegrationPoints(); it++){
        
       xsi[0] = nQuad.PointList(index,0);
       xsi[1] = nQuad.PointList(index,1);
            
    //    //Computes the velocity shape functions
    //    shapeQuad.Shape(xsi,phi_);

       for (int j=0; j<fMesh->NElNodes(); j++){
           fIntPointDistFunction[index] +=  fMesh->getNumericalIntegration()->phi_(j,index) * nodalval[j];
       };
       // intPointWeightFunction(index) = 1.;
       index++;
    }; 

    index = 0;
     return;
};

//------------------------------------------------------------------------------
//------------------COMPUTES THE INTEGRATION POINT COORDINATE-------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::getIntegPointCoordinates(){
    int DIM = tshape::Dimension;
    DEG = fMesh->GetDefaultOrder();

    IntegQuadrature sQuad(DIM,DEG);
    fIntPointCoordinates.resize(sQuad.getNumberOfIntegrationPoints(),2);
    fIntPointCoordinates.setZero();
    
    VecDouble xsi(DIM);
    ShapeFunction shapeQuad(DIM,DEG);
    VecDouble phi_(fMesh->NElNodes());
    fIntPointCoordinates.resize(sQuad.getNumberOfIntegrationPoints(),DIM);

    for (int i = 0; i < sQuad.getNumberOfIntegrationPoints(); i++){
        double x[DIM] = {};

        for (int k = DIM; k--; ) xsi[k] = sQuad.PointList(i,k);

        shapeQuad.Shape(xsi,phi_);

        for (int k = DIM; k--; ) fIntPointCoordinates(i,k) = 0.;

        for (int j = 0; j < fMesh->NElNodes(); j++)
            for (int k = DIM; k--; )
                fIntPointCoordinates(i,k) += fMesh->NodeVec()[fConnect[j]] -> getCoordinateValue(k) * phi_[j];
        
    };

    return;
};



//------------------------------------------------------------------------------
//---------------------------CLEAR ELEMENT VARIABLES----------------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::clearVariables(){

    int DIM = tshape::Dimension;
    IntegQuadratureSpecial sQuad(DIM,DEG);

    for (int i=0; i < sQuad.getNumberOfIntegrationPoints(); i++){
        intPointWeightFunction[i] = 1.;
        for (int j =0; j<DIM; j++) {
            fIntPointCoordinates(i,j) = 0.;
        }
    };
    
    getIntegPointCoordinates();

    return;
}; 

//------------------------------------------------------------------------------
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::setIntersectionParameters(VecDouble &x, VecDouble &X) {
    xK.resize(2);
    XK.resize(2);
    xK[0] = x[0]; xK[1] = x[1]; 
    XK[0] = X[0]; XK[1] = X[1];
    return;
};


//------------------------------------------------------------------------------
//-------------------------SPATIAL TRANSFORM - JACOBIAN-------------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::ComputeJacobian(int index) {

    int DIM = tshape::Dimension;
    fIntegData.fA0Inv.resize(DIM,DIM);
    fIntegData.fA0.resize(DIM,DIM);
    fIntegData.fA0.setZero();
    fIntegData.fX.resize(DIM);
    fIntegData.fX.setZero();

    double &alpha_f = fMesh->getProblemParameters().getAlphaF();

    fIntegData.fA0.setZero();
    for (int i = fMesh->NElNodes(); i--; ){
        for (int j = DIM; j--; ){
            // Approximate the integration space
            fIntegData.fX[j] = fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j) ;
            // xna_[j] = alpha_f * fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j) + 
            //           (1. - alpha_f) * fMesh->NodeVec()[fConnect[i]] -> getPreviousCoordinateValue(j);
            for (int k = DIM; k--; ){
                fIntegData.fA0(j,k) += fIntegData.fX[j] * fMesh->getNumericalIntegration()->dphi_[i](k,index);
                // dx_dxsi(j,k) += xna_[j] * dphi(i,k);
            };
        };
    };

    //Computing the jacobian determinant and Inverse
    fIntegData.fJacA0 = fIntegData.fA0.determinant();
    fIntegData.fA0Inv = fIntegData.fA0.inverse().transpose();

    return;
};

template<class tshape>
void ElementT<tshape>::ComputeCurrentJacobian(int index) {

    int DIM = tshape::Dimension;
    fIntegData.fA1.resize(DIM,DIM);
    fIntegData.fA1.setZero();
    fIntegData.fX.resize(DIM);
    fIntegData.fX.setZero();

    double &alpha_f = fMesh->getProblemParameters().getAlphaF();

    fIntegData.fA1.setZero();
    for (int i = fMesh->NElNodes(); i--; ){
        for (int j = DIM; j--; ){
            // Approximate the integration space
            fIntegData.fX[j] = fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j) + fMesh->NodeVec()[fConnect[i]] -> GetSolution(j) ;
            // xna_[j] = alpha_f * fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j) + 
            //           (1. - alpha_f) * fMesh->NodeVec()[fConnect[i]] -> getPreviousCoordinateValue(j);
            for (int k = DIM; k--; ){
                fIntegData.fA1(j,k) += fIntegData.fX[j] * fMesh->getNumericalIntegration()->dphi_[i](k,index);
                // dx_dxsi(j,k) += xna_[j] * dphi(i,k);
            };
        };
    };

    //Computing the jacobian determinant and Inverse
    fIntegData.fJacA1 = fIntegData.fA1.determinant();

    return;
};

//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::ComputeSpatialDerivatives() {
    
    // typename QuadShapeFunction<2,2>::ValueDDeriv ddphi;
    int DIM = tshape::Dimension;
    MatrixDouble dphi(fMesh->NElNodes(),DIM);
    
    ShapeFunction shapeQuad(DIM,DEG);
    
    shapeQuad.ShapeGradient(fIntegData.fAdimCoord,dphi);
    // shapeQuad.ShapeHessian(xsi,ddphi);
    
    fIntegData.fDPhiX0.setZero();

    //Shape functions spatial first derivatives
    fIntegData.fDPhiX0 = dphi * fIntegData.fA0Inv.transpose();

    return;
};

template<class tshape>
void ElementT<tshape>::ComputeCurrentSpatialDerivatives() {
    
    // typename QuadShapeFunction<2,2>::ValueDDeriv ddphi;
    int DIM = tshape::Dimension;
    MatrixDouble dphi(fMesh->NElNodes(),DIM);
    
    ShapeFunction shapeQuad(DIM,DEG);
    
    shapeQuad.ShapeGradient(fIntegData.fAdimCoord,dphi);
    // shapeQuad.ShapeHessian(xsi,ddphi);
    
    fIntegData.fDPhiX1.setZero();

    //Shape functions spatial first derivatives
    fIntegData.fDPhiX1 = dphi * fIntegData.fA1.inverse().transpose();

    return;
};
//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::getHighOrderSpatialDerivatives(VecDouble &xsi, MatrixDouble &ainv_, MatrixDouble &dphi_dx, MatrixDouble &dDphi_dx) {
    
    int DIM = tshape::Dimension;
    dDphi_dx.setZero();
    std::vector<MatrixDouble> ddphi(fMesh->NElNodes(),MatrixDouble(DIM,DIM));
    for (int i = 0; i < fMesh->NElNodes(); i++)
    {
        ddphi[i].setZero();
    }
    

    ShapeFunction shapeQuad(DIM,DEG);
    shapeQuad.ShapeHessian(xsi,ddphi);
    
    //These derivatives are computed with basis in this reference:
    //https://scicomp.stackexchange.com/questions/25196/implementing-higher-order-derivatives-for-finite-element
    MatrixDouble matAux,invSecDeriv;
    VecDouble HODerivatives, vecAux;
    // if (!isSecondDerivativeInverted){
        if (DIM == 2){
            matAux.resize(3,3); HODerivatives.resize(3); vecAux.resize(3); invSecDeriv.resize(3,3);
            matAux.setZero(); HODerivatives.setZero(); vecAux.setZero(); invSecDeriv.setZero();

            matAux(0,0) = ainv_(0,0) * ainv_(0,0);
            matAux(0,1) = ainv_(0,1) * ainv_(0,1);//Pode estar errado, e ser(1,0)ou seja a transposta
            matAux(0,2) = 2. * ainv_(0,0) * ainv_(0,1);
            
            matAux(1,0) = ainv_(1,0) * ainv_(1,0);
            matAux(1,1) = ainv_(1,1) * ainv_(1,1);
            matAux(1,2) = 2. * ainv_(1,0) * ainv_(1,1);
            
            matAux(2,0) = ainv_(0,0) * ainv_(1,0);
            matAux(2,1) = ainv_(0,1) * ainv_(1,1);
            matAux(2,2) = ainv_(0,0) * ainv_(1,1) + ainv_(1,0) * ainv_(0,1);

            invSecDeriv = matAux.inverse();

        } else if (DIM == 3){
            matAux.resize(6,6); HODerivatives.resize(6); vecAux.resize(6);invSecDeriv.resize(6,6);
            matAux.setZero(); HODerivatives.setZero(); vecAux.setZero(); invSecDeriv.setZero();
        } else {
            PanicButton();
        }
    // }
    
    //Shape functions spatial second derivatives
    double ddx_dxsi, ddx_deta, ddx_dxsideta, ddy_dxsi, ddy_deta, ddy_dxsideta;
    double &alpha_f = fMesh->getProblemParameters().getAlphaF();
    VecDouble xna_(DIM);

    for (int i = fMesh->NElNodes(); i--; ){
        xna_.setZero();
        for (int j = DIM; j--; ){
            // Approximate the integration space
            xna_[j] = alpha_f * fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j) + 
                      (1. - alpha_f) * fMesh->NodeVec()[fConnect[i]] -> getPreviousCoordinateValue(j);
        }
        ddx_dxsi += xna_[0] * ddphi[i](0,0);
        ddx_deta += xna_[0] * ddphi[i](1,1);
        ddx_dxsideta += xna_[0] * ddphi[i](0,1);
        ddy_dxsi += xna_[1] * ddphi[i](0,0);
        ddy_deta += xna_[1] * ddphi[i](1,1);
        ddy_dxsideta += xna_[1] * ddphi[i](0,1);
    };

    if (DIM == 2){
        VecDouble vecAux2(3);
        for (int i = fMesh->NElNodes(); i--; ){
            vecAux[0] = ddphi[i](0,0) - dphi_dx(i,0)*ddx_dxsi - dphi_dx(i,1)*ddy_dxsi;
            vecAux[1] = ddphi[i](1,1) - dphi_dx(i,0)*ddx_deta - dphi_dx(i,1)*ddy_deta;
            vecAux[2] = ddphi[i](0,1) - dphi_dx(i,0)*ddx_dxsideta - dphi_dx(i,1)*ddy_dxsideta;

            vecAux2 = invSecDeriv * vecAux;

            dDphi_dx(i,0) = vecAux2[0];
            dDphi_dx(i,1) = vecAux2[1];
            dDphi_dx(i,2) = vecAux2[2];
        }
        
        
    }

    return;
};


//------------------------------------------------------------------------------
//---------------INTERPOLATES MESH VELOCITY AND ITS DERIVATIVES-----------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::interpolateMeshVelocity(int &index, VecDouble &umesh_, VecDouble &umeshPrev_) {

    int DIM = tshape::Dimension;
    umesh_.setZero();
    umeshPrev_.setZero();
    
    for (int i = fMesh->NElNodes(); i--; ){
        double shapeFi = fMesh->getNumericalIntegration()-> phi_(i,index);
        for (int j = DIM; j--; ){
            umesh_[j] += fMesh->NodeVec()[fConnect[i]] -> getMeshVelocity(j) * shapeFi;
        }
    }

    return;
}

template<class tshape>
void ElementT<tshape>::interpolateSolution(int &index, VecDouble &u_) {
    u_.setZero();
    for (int i = fMesh->NElNodes(); i--; ){
        double shapeFi = fMesh->getNumericalIntegration()-> phi_(i,index);
        int nstate = fMesh->NodeVec()[fConnect[i]]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            u_[j] += fMesh->NodeVec()[fConnect[i]] -> GetSolution(j) * shapeFi;
        }
    }
}

template<class tshape>
void ElementT<tshape>::interpolateSolution(VecDouble &phi, VecDouble &u_) {
    u_.setZero();
    for (int i = fMesh->NElNodes(); i--; ){
        double shapeFi = phi(i);
        int nstate = fMesh->NodeVec()[fConnect[i]]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            u_[j] += fMesh->NodeVec()[fConnect[i]] -> GetSolution(j) * shapeFi;
        }
    }
}

//------------------------------------------------------------------------------
//----------------------------INTERPOLATES VELOCITY-----------------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::interpolateSolDerivatives(MatrixDouble &du_dx) {
    du_dx.setZero();    
    int DIM = tshape::Dimension;
    for (int i = fMesh->NElNodes(); i--; ){
        int nstate = fMesh->NodeVec()[fConnect[i]]->GetNStateVariables();
        for (int j = DIM; j--; ){
            for (int k = nstate; k--; ){
                du_dx(k,j) += fMesh->NodeVec()[fConnect[i]] -> GetSolution(k) * fIntegData.fDPhiX0(i,j);
            }
        }
    }
}

//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::getBoundaryLoad(VecDouble &xsi, VecDouble &load) {
    // std::cout << "asdasd 0 " << std::endl;
    // int nBdNodes = fMesh->NBdNodes() = 3*(1-DEG)+DIM*(2*DEG-1);

    // double &visc_ = fMesh->getProblemParameters().GetViscosity();
    // double &alpha_f = fMesh->getProblemParameters().getAlphaF();

    // VecDouble phi_(fMesh->NElNodes());

    // MatrixDouble dphi_dx(fMesh->NElNodes(),DIM);

    // double shearStress[DIM][DIM] = {};
    
    // MatrixDouble dphi(fMesh->NElNodes(),DIM);
    // MatrixDouble ainv_(DIM,DIM);

    // double ident[DIM][DIM] = {}; ident[0][0] = 1.; ident[1][1] = 1.;
    // ShapeFunction shapeQuad(DIM,DEG);
    
    // // std::cout << "asdasd 1 " << std::endl;
    // fIntegData.fAdimCoord = xsi;
    // //Computes the shape functions        
    // shapeQuad.Shape(xsi,phi_);
    // double djac_ = 0.;
    // //Computes the jacobian matrix
    // ComputeJacobian(0);
    // //Computes spatial derivatives
    // ComputeSpatialDerivatives();

    // shapeQuad.ShapeGradient(xsi, dphi); 

    // //Velocity Derivatives
    // MatrixDouble du_dx(DIM,DIM), duprev_dx(DIM,DIM), duna_dx(DIM,DIM);
    // // interpolateVelDerivatives(dphi_dx, du_dx, duprev_dx);
    // duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
    //     // std::cout << "asdasd 4 " << std::endl;
    // //Pressure
    // double p_;
    // VecDouble dp_dx(DIM);
    // int index = 0;
    // // interpolatePressure(index, dphi_dx, p_, dp_dx);
    // // std::cout << "asdasd 5 " << std::endl;
    // shearStress[0][0] = 2. * visc_ * duna_dx(0,0);
    // shearStress[0][1] = visc_ * (duna_dx(0,1) + duna_dx(1,0));
    // shearStress[1][0] = visc_ * (duna_dx(0,1) + duna_dx(1,0));
    // shearStress[1][1] = 2. * visc_ * duna_dx(1,1);

    // int nodesb_[fMesh->NBdNodes()];
    // VecDouble xsiB(DIM-1);
    
    // if (DEG == 2){
    //     if(fSideInBoundary == 0){
    //         nodesb_[0] = fConnect[1]; 
    //         nodesb_[1] = fConnect[4]; 
    //         nodesb_[2] = fConnect[2];
    //         xsiB[0] = 2. * xsi[0] - 1.;
    //     }else{
    //         if(fSideInBoundary == 1){
    //             nodesb_[0] = fConnect[2]; 
    //             nodesb_[1] = fConnect[5]; 
    //             nodesb_[2] = fConnect[0]; 
    //             xsiB[0] = 2. * xsi[1] - 1.;
    //         }else{
    //             nodesb_[0] = fConnect[0];
    //             nodesb_[1] = fConnect[3];
    //             nodesb_[2] = fConnect[1];
    //             xsiB[0] = 1. - 2. * xsi[0];
    //         };        
    //     };
    // } else {
    //     if(fSideInBoundary == 0){
    //         nodesb_[0] = fConnect[2]; 
    //         nodesb_[1] = fConnect[1];
    //         xsiB[0] = 2. * xsi[0] - 1.;
    //     }else{
    //         if(fSideInBoundary == 1){
    //             nodesb_[0] = fConnect[0]; 
    //             nodesb_[1] = fConnect[2];
    //             xsiB[0] = 2. * xsi[1] - 1.;
    //         }else{
    //             nodesb_[0] = fConnect[1];
    //             nodesb_[1] = fConnect[0];
    //             xsiB[0] = 1. - 2. * xsi[0];
    //         };        
    //     };
    // }

    // VecDouble phib_(fMesh->NBdNodes());

    // MatrixDouble dphib_(fMesh->NBdNodes(),DIM-1);

    // double dx_dxsiB[3][DIM-1] = {};
    // double xna_[3] = {};

    // BoundShapeFunction::getShapeFunction(DIM,DEG,xsiB,phib_,dphib_);

    // // shapeBound.getShapeFunction(xsiB,phib_,dphib_);

    // // std::cout << "asdasd 7 " << std::endl;
    // for (int i = 0; i < fMesh->NBdNodes(); i++){
    //     for (int j = 0; j < DIM; j++){
    //         // Approximate the integration space
    //         xna_[j] = alpha_f * fMesh->NodeVec()[nodesb_[i]] -> getCoordinateValue(j) + 
    //                 (1.-alpha_f) * fMesh->NodeVec()[nodesb_[i]] -> getPreviousCoordinateValue(j);
    //         for (int k = DIM-1; k--; ) dx_dxsiB[j][k] += xna_[j] * dphib_(i,k);
    //     };
    // };
    // // std::cout << "asdasd 8 " << std::endl;
    // double Maux[DIM-1][DIM-1] = {};
    //     for (int i = DIM-1; i--; )
    //         for (int j = DIM-1; j--; )
    //             for (int k = 3; k--; )
    //                 Maux[i][j] += dx_dxsiB[k][i] * dx_dxsiB[k][j];

    // double djacb_ = 0.;
    // if (DIM == 2){
    //     djacb_ = std::sqrt(Maux[0][0]);
    // } else {
    //     djacb_ = std::sqrt(Maux[0][0] * Maux[1][1] - Maux[0][1] * Maux[1][0]);
    // }

    // VecDouble n_vector(DIM);
    // if (DIM == 2){
    //     n_vector[0] =  dx_dxsiB[1][0] / djacb_;
    //     n_vector[1] = -dx_dxsiB[0][0] / djacb_;
    // } else {
    //     n_vector[0] = (dx_dxsiB[1][0]*dx_dxsiB[2][1] - dx_dxsiB[2][0]*dx_dxsiB[1][1]) / djacb_;
    //     n_vector[1] = (dx_dxsiB[2][0]*dx_dxsiB[0][1] - dx_dxsiB[0][0]*dx_dxsiB[2][1]) / djacb_;
    //     n_vector[2] = (dx_dxsiB[0][0]*dx_dxsiB[1][1] - dx_dxsiB[1][0]*dx_dxsiB[0][1]) / djacb_;
    // }

    // // for (int i = 0; i < fMesh->NBdNodes(); i++) fMesh->NodeVec()[nodesb_[i]] -> setInnerNormal(n_vector);
    
    // // std::cout << "asdasd 9 " << std::endl;
    // // std::cout << "dx_dxsiB " << dx_dxsiB[0][0] << " " << dx_dxsiB[1][0] << " " << dx_dxsiB[2][0] << std::endl;
    // // std::cout << "N vector " << xna_[0] << " " << xna_[1] << " " << n_vector[0] << " " << n_vector[1] << " " << djacb_<< " " << fIndex << std::endl;

    // // if (fSideInBoundary == 1){
    // //     for (int i = 0; i < 6; i++){
    // //         t_vector[0] -= dphi[1][i] * fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(0);
    // //         t_vector[1] -= dphi[1][i] * fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(1);
    // //     };        
    // // };

    // // if (fSideInBoundary == 2){
    // //     for (int i = 0; i < 6; i++){
    // //         t_vector[0] += dphi[0][i] * fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(0);
    // //         t_vector[1] += dphi[0][i] * fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(1);
    // //     };        
    // // };

    // // if (fSideInBoundary == 0){
    // //     std::cout << "VERIFICAR VETOR NORMAL - getBoundaryLoad" << std::endl;
    // // };

    // // n_vector[0] =  t_vector[1] / std::sqrt(t_vector[0]*t_vector[0] + t_vector[1]*t_vector[1]);
    // // n_vector[1] = -t_vector[0] / std::sqrt(t_vector[0]*t_vector[0] + t_vector[1]*t_vector[1]);

    // for (int i = 0; i < DIM; i++)
    //     for (int j = 0; j < DIM; j++)
    //         load[i] -= -p_ * ident[i][j] * n_vector[j] + shearStress[i][j] * n_vector[j];
    // // load = -p_ * prod(ident,n_vector) + prod(shearStress,n_vector);

    // //std::cout << "N Vector " << fSideInBoundary << " " <<  n_vector(0) << " " << n_vector(1) << " " << load(0) << " " << load(1) << " " << p_ << std::endl;
    
    // // std::cout << "asdasd 10 " << std::endl;


    return;
};


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//----------------MOUNTS EACH TYPE OF INCOMPRESSIBLE FLOW PROBEM----------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

    int DIM = tshape::Dimension;
    DEG = fMesh->GetDefaultOrder();
    
    ShapeFunction shapeQuad(DIM,DEG);
    int index = 0;
    IntegQuadrature nQuad(DIM,DEG);
    fIntegData.fAdimCoord.resize(DIM);
    for(int it = 0; it < nQuad.getNumberOfIntegrationPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) fIntegData.fAdimCoord[k] = nQuad.PointList(index,k);

        //Returns the quadrature integration weight
        fIntegData.fWeight = nQuad.WeightList(index);

        //Computes the jacobian matrix
        ComputeJacobian(index);

        //Computes spatial derivatives
        ComputeSpatialDerivatives();

        //Computes the element diffusion/viscosity matrix
        ComputeStiffness(index, jacobianNRMatrix);

        //Computes the RHS vector
        ComputeResidual(index, rhsVector); 

        index++;        
    };  
    // std::cout << "\nStiffness Element " << this->Index() << "\n" << jacobianNRMatrix;
    // std::cout << "\nrhsVector Element " << this->Index() << "\n" << rhsVector;
    //Apply boundary conditions
    ApplyBC(jacobianNRMatrix, rhsVector);

    return;
};


//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::ComputeElContribution(std::vector<MatrixDouble> &jacobianNRMatrix, std::vector<VecDouble> &rhsVector){

    int DIM = tshape::Dimension;
    DEG = fMesh->GetDefaultOrder();
    fIntegData.fA0Inv.resize(DIM,DIM);
    fIntegData.fAdimCoord.resize(DIM);

    ShapeFunction shapeQuad(DIM,DEG);
    int index = 0;
    IntegQuadrature nQuad(DIM,DEG);

    for(int it = 0; it < nQuad.getNumberOfIntegrationPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) fIntegData.fAdimCoord[k] = nQuad.PointList(index,k);

        //Returns the quadrature integration weight
        fIntegData.fWeight = nQuad.WeightList(index);

        //Computes the jacobian matrix
        ComputeJacobian(index);

        //Computes spatial derivatives
        ComputeSpatialDerivatives();

        //Computes the element diffusion/viscosity matrix
        ComputeStiffness(index, jacobianNRMatrix);

        //Computes the RHS vector
        ComputeResidual(index, rhsVector); 

        index++;        
    };  
    // std::cout << "\nStiffness Element " << this->Index() << "\n" << jacobianNRMatrix[1];
    // std::cout << "\nrhsVector Element " << this->Index() << "\n" << rhsVector[1];
    //Apply boundary conditions
    ApplyBC(jacobianNRMatrix, rhsVector);

    return;
};

#include "ShapeHexahedron.h"
#include "ShapeOneD.h"
#include "ShapeQuadrilateral.h"
#include "ShapePoint.h"
#include "ShapeTetrahedron.h"
#include "ShapeTriangle.h"

template class ElementT<ShapePoint>;
template class ElementT<ShapeOneD>;
template class ElementT<ShapeTriangle>;
template class ElementT<ShapeQuadrilateral>;
template class ElementT<ShapeTetrahedron>;
template class ElementT<ShapeHexahedron>;