#include "ElementT.h"
#include "PositionalTruss.h"
#include "ElasticityPositional2D.h"
#include "L2Projection.h"
#include "PositionalFrame2D.h"

template<class compshape>
ElementT<compshape>::ElementT() : Element(){
};

template<class compshape>
ElementT<compshape>::ElementT(int64_t index, GeoElement* gel, CompMesh* mesh, WeakForm *wf) : Element(){
    fReference = gel;
    int DIM = compshape::Dimension;
    fMesh = mesh;
    fConnect.resize(compshape::NSides);
    fIndex = index;
    fWeakForm = wf;
    if (fWeakForm) nLocDOF = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()) * fWeakForm->NState();

    
    int increase = 0;
    if(wf->GetExactSolution()) increase = 2;

    fIntRule.SetOrder(2*this->fMesh->GetDefaultOrder()+increase);
    // fIntRule.SetOrder(2);

    fIntegData.fWeightFunction.resize(fIntRule.NPoints());
    fIntegData.fDistFunction.resize(fIntRule.NPoints());
    fIntegData.fPrevWeightFunction.resize(fIntRule.NPoints());

    fIntegData.fWeightFunction.fill(1.);
    fIntegData.fPrevWeightFunction.fill(1.);

    fIntegData.fAdimCoord.resize(DIM);

    fIntegData.fA0Inv.resize(DIM,DIM);
    fIntegData.fA0Inv.setZero();
    fIntegData.fAxes0.resize(3,DIM);
    fIntegData.fAxes0.setZero();
    fIntegData.fA0.resize(DIM,DIM);
    fIntegData.fA0.setZero();
    fIntegData.fX.resize(3);
    fIntegData.fX.setZero();

    fIntegData.fPhi.resize(compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()));
    fIntegData.fPhi.setZero();

    fIntegData.fDPhi.resize(compshape::Dimension,compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()));
    fIntegData.fDPhi.setZero();

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template<class compshape>
void ElementT<compshape>::ComputeSpatialDerivatives() {
    fIntegData.fPhi.setZero();
    fIntegData.fDPhi.setZero();

    VecInt orders(compshape::NSides);
    if (fConnect.size() != compshape::NSides) {
        orders.setZero();
    } else {
        for (int i = compshape::NSides; i--; ) orders[i] = fConnect[i] -> GetOrder();
    }

    //Shape functions
    compshape::Shape(fIntegData.fAdimCoord,fIntegData.fPhi,fIntegData.fDPhi,orders);
    // shapeQuad.ShapeHessian(xsi,ddphi);

    //Shape functions spatial first derivatives
    fIntegData.fDPhiX0 = fIntegData.fA0Inv.transpose() * fIntegData.fDPhi;

    return;
};

template<class compshape>
void ElementT<compshape>::ComputeCurrentSpatialDerivatives() {
    VecInt orders(compshape::NSides);
    for (int i = compshape::NSides; i--; ) orders[i] = fConnect[i] -> GetOrder();

    MatrixDouble DphiComp(compshape::Dimension,compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()));
    compshape::Shape(fIntegData.fAdimCoord,fIntegData.fPhi,DphiComp,orders);
    // shapeQuad.ShapeHessian(xsi,ddphi);
    
    fIntegData.fDPhiX1.setZero();

    //Shape functions spatial first derivatives
    fIntegData.fDPhiX1 = fIntegData.fA1.inverse().transpose() * DphiComp;

    return;
};
//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template<class compshape>
void ElementT<compshape>::ComputeHighOrderSpatialDerivatives() {
    
    // int DIM = compshape::Dimension;
        
    // int nshape = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder());

    // dDphi_dx.setZero();
    // std::vector<MatrixDouble> ddphi(nshape,MatrixDouble(DIM,DIM));
    // for (int i = 0; i < nshape; i++)
    // {
    //     ddphi[i].setZero();
    // }
    
    // compshape::ShapeHessian(xsi,ddphi);
    
    // //These derivatives are computed with basis in this reference:
    // //https://scicomp.stackexchange.com/questions/25196/implementing-higher-order-derivatives-for-finite-element
    // MatrixDouble matAux,invSecDeriv;
    // VecDouble HODerivatives, vecAux;
    // // if (!isSecondDerivativeInverted){
    //     if (DIM == 2){
    //         matAux.resize(3,3); HODerivatives.resize(3); vecAux.resize(3); invSecDeriv.resize(3,3);
    //         matAux.setZero(); HODerivatives.setZero(); vecAux.setZero(); invSecDeriv.setZero();

    //         matAux(0,0) = ainv_(0,0) * ainv_(0,0);
    //         matAux(0,1) = ainv_(0,1) * ainv_(0,1);//Pode estar errado, e ser(1,0)ou seja a transposta
    //         matAux(0,2) = 2. * ainv_(0,0) * ainv_(0,1);
            
    //         matAux(1,0) = ainv_(1,0) * ainv_(1,0);
    //         matAux(1,1) = ainv_(1,1) * ainv_(1,1);
    //         matAux(1,2) = 2. * ainv_(1,0) * ainv_(1,1);
            
    //         matAux(2,0) = ainv_(0,0) * ainv_(1,0);
    //         matAux(2,1) = ainv_(0,1) * ainv_(1,1);
    //         matAux(2,2) = ainv_(0,0) * ainv_(1,1) + ainv_(1,0) * ainv_(0,1);

    //         invSecDeriv = matAux.inverse();

    //     } else if (DIM == 3){
    //         matAux.resize(6,6); HODerivatives.resize(6); vecAux.resize(6);invSecDeriv.resize(6,6);
    //         matAux.setZero(); HODerivatives.setZero(); vecAux.setZero(); invSecDeriv.setZero();
    //     } else {
    //         PanicButton();
    //     }
    // // }
    
    // //Shape functions spatial second derivatives
    // double ddx_dxsi, ddx_deta, ddx_dxsideta, ddy_dxsi, ddy_deta, ddy_dxsideta;
    // VecDouble xna_(DIM);

    // for (int i = nshape; i--; ){
    //     xna_.setZero();
    //     for (int j = DIM; j--; ){
    //         // Approximate the integration space
    //         xna_[j] = fMesh->Reference()->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j);
    //     }
    //     ddx_dxsi += xna_[0] * ddphi[i](0,0);
    //     ddx_deta += xna_[0] * ddphi[i](1,1);
    //     ddx_dxsideta += xna_[0] * ddphi[i](0,1);
    //     ddy_dxsi += xna_[1] * ddphi[i](0,0);
    //     ddy_deta += xna_[1] * ddphi[i](1,1);
    //     ddy_dxsideta += xna_[1] * ddphi[i](0,1);
    // };

    // if (DIM == 2){
    //     VecDouble vecAux2(3);
    //     for (int i = nshape; i--; ){
    //         vecAux[0] = ddphi[i](0,0) - dphi_dx(i,0)*ddx_dxsi - dphi_dx(i,1)*ddy_dxsi;
    //         vecAux[1] = ddphi[i](1,1) - dphi_dx(i,0)*ddx_deta - dphi_dx(i,1)*ddy_deta;
    //         vecAux[2] = ddphi[i](0,1) - dphi_dx(i,0)*ddx_dxsideta - dphi_dx(i,1)*ddy_dxsideta;

    //         vecAux2 = invSecDeriv * vecAux;

    //         dDphi_dx(i,0) = vecAux2[0];
    //         dDphi_dx(i,1) = vecAux2[1];
    //         dDphi_dx(i,2) = vecAux2[2];
    //     }
        
        
    // }

    return;
};


template<class compshape>
void ElementT<compshape>::setIntegPointWeightFunction() {
    //It produces wrong results for constant weight functions
    VecDouble xsi(Dimension());    

    for(int i = 0; i < fIntRule.NPoints(); i++) {
        fIntegData.fPrevWeightFunction[i] = fIntegData.fWeightFunction[i];
        fIntegData.fWeightFunction[i] = 0.;
    }

    int index=0;

    for(int it = 0; it < fIntRule.NPoints(); it++){
        
       xsi[0] = fIntRule.PointList(index,0);
       xsi[1] = fIntRule.PointList(index,1);

       for (int j=0; j<compshape::NSides; j++){
        std::cout << "alterar armazenamento do weightfunction para vector de connects" << std::endl;
        PanicButton();
           fIntegData.fWeightFunction[index] += fIntegData.fPhi[j] * fConnect[j] -> getWeightFunction();
       };
       // fIntegData.fWeightFunction(index) = 1.;
       index++;
    }; 

    // index = 0;
     return;
};

// //------------------------------------------------------------------------------
// //----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
// //------------------------------------------------------------------------------


template<class compshape>
void ElementT<compshape>::ComputeIntPointDistFunction(VecDouble &nodalval) {
    
    PanicButton();
    // int DIM = Dimension();
    // fIntegData.fAdimCoord.resize(DIM);
    // fIntegData.fPhi.resize(geoshape::NShape);
    
    // // for(int i = 0; i < nQuad.getNumberOfIntegrationPoints(); i++) {
    // //     intPointWeightFunctionPrev[i] = intPointWeightFunction[i];
    // //     intPointWeightFunction[i] = 0.;
    // // }

    // int index=0;
    // fIntegData.fDistFunction.setZero();

    // VecInt orders(compshape::NSides);
    // for (int i = compshape::NSides; i--; ) orders[i] = this->fMesh->ConnectVec()[fConnect[i]] -> GetOrder();

    // for(int it = 0; it < fIntRule.NPoints(); it++){
        
    //     for (int i=0; i<DIM; i++) fIntegData.fAdimCoord[i] = fIntRule.PointList(index,i);
           
    //     compshape::Shape(fIntegData.fAdimCoord,fIntegData.fPhi,fIntegData.fDPhi,orders);

    //     for (int j=0; j<geoshape::NShape; j++){
    //             fIntegData.fDistFunction[index] +=  fIntegData.fPhi[j] * nodalval[j];
    //     };

    //     // intPointWeightFunction(index) = 1.;
    //     index++;
    // }; 

    // index = 0;
    return;
};

//------------------------------------------------------------------------------
//------------------COMPUTES THE INTEGRATION POINT COORDINATE-------------------
//------------------------------------------------------------------------------
template<class compshape>
void ElementT<compshape>::ComputeIntegPointCoordinates(){

    int DIM = compshape::Dimension;
    fIntPointCoordinates.resize(fIntRule.NPoints(),3);
    fIntPointCoordinates.setZero();

    int nshape = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder());
    VecDouble xsi(DIM);
    VecDouble phi_(nshape);
    MatrixDouble dphi_(DIM,nshape);

    VecInt orders(compshape::NSides);
    for (int i = compshape::NSides; i--; ) orders[i] = fConnect[i] -> GetOrder();

    VecInt fGeoNodes = fReference->getGeometricNodes();

    for (int i = 0; i < fIntRule.NPoints(); i++){
        for (int k = DIM; k--; ) xsi[k] = fIntRule.PointList(i,k);

        compshape::Shape(xsi,phi_,dphi_,orders);

        for (int j = 0; j < nshape; j++)
            for (int k = 0; k < 3; k++)
                fIntPointCoordinates(i,k) += fMesh->Reference()->NodeVec()[fGeoNodes[j]] -> getCoordinateValue(k) * phi_[j];
        
    };

    return;
};





//------------------------------------------------------------------------------
//---------------INTERPOLATES MESH VELOCITY AND ITS DERIVATIVES-----------------
//------------------------------------------------------------------------------
template<class compshape>
void ElementT<compshape>::interpolateMeshVelocity(int &index, VecDouble &umesh_, VecDouble &umeshPrev_) {

    int DIM = compshape::Dimension;
    umesh_.setZero();
    umeshPrev_.setZero();

    PanicButton();
    
    // for (int i = compshape::NShape; i--; ){
    //     double shapeFi = fIntegData.fPhi[i];
    //     for (int j = DIM; j--; ){
    //         umesh_[j] += fMesh->ConnectVec()[fConnect[i]] -> getMeshVelocity(j) * shapeFi;
    //     }
    // }

    return;
}

template<class compshape>
void ElementT<compshape>::interpolateSolution(int &index, VecDouble &u_) {
    u_.setZero();

    for (int i = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()); i--; ){
        double shapeFi = fIntegData.fPhi[i];
        int nstate = fConnect[i]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            u_[j] += fConnect[i] -> GetSolution(j) * shapeFi;
        }
    }
}
template<class compshape>
double ElementT<compshape>::InterpolateVariable(VecDouble &nValues, int point) {
    double val = 0.;
    int nshape = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder());
    fIntegData.fPhi.resize(nshape);
    fIntegData.fPhi.setZero();
    fIntegData.fAdimCoord.resize(compshape::Dimension);
    for (int i=0; i<compshape::Dimension; i++) fIntegData.fAdimCoord[i] = fIntRule.PointList(point,i);

    VecInt orders(compshape::NSides);
    for (int i = compshape::NSides; i--; ) orders[i] = fConnect[i] -> GetOrder();

    compshape::Shape(fIntegData.fAdimCoord,fIntegData.fPhi,fIntegData.fDPhi,orders);
    for (int i = nshape; i--; ){
        double shapeFi = fIntegData.fPhi[i];
        val += nValues[i] * shapeFi;
    }
    return val;
}

template<class compshape>
void ElementT<compshape>::interpolateSolDTimeDerivatives() {
    fIntegData.fDSolDt.setZero();
    fIntegData.fDSolDDt.setZero();

    int count = 0;
    for (int iside = 0; iside < compshape::NSides; iside++){
        
        int nshape = compshape::NShapeFunctions(iside,this->fMesh->GetDefaultOrder());
        if (nshape == 0) continue;
        int nstate = fMesh->NState();
    
        for (int ishape = 0; ishape < nshape; ishape++){
            double shapeFi = fIntegData.fPhi[count];
            for (int j = 0; j < nstate; j++ ){
                fIntegData.fDSolDt[j] += fConnect[iside] -> GetDSolutionDTime(j) * shapeFi;
                fIntegData.fDSolDDt[j] += fConnect[iside] -> GetDSolutionDDTime(j) * shapeFi;
            }
            count++;
        }
    }
}

template<class compshape>
void ElementT<compshape>::interpolateSolDTimeDerivatives(VecDouble &du_dt, VecDouble &du_ddt) {
    fIntegData.fDSolDt.setZero();
    fIntegData.fDSolDDt.setZero();
    for (int i = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()); i--; ){
        double shapeFi = fIntegData.fPhi[i];
        int nstate = fConnect[i]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            fIntegData.fDSolDt[j] += fConnect[i] -> GetDSolutionDTime(j) * shapeFi;
            fIntegData.fDSolDDt[j] += fConnect[i] -> GetDSolutionDDTime(j) * shapeFi;
        }
    }
    du_dt = fIntegData.fDSolDt;
    du_ddt = fIntegData.fDSolDDt;
}


template<class compshape>
void ElementT<compshape>::interpolateSolution() {
    fIntegData.fSol.setZero();
    int count = 0;
    if (fIntegData.fSolPrev.size() != 0) fIntegData.fSolPrev.setZero();
    if (fConnect.size() == 0) return;

    if (fConnect.size()!=compshape::NSides) {
        int nstate = fMesh->NState();
        for (int j = 0; j < nstate; j++ ){
            fIntegData.fSol[j] += fConnect[0] -> GetSolution(j);
            if (fIntegData.fSolPrev.size() != 0) {
                fIntegData.fSolPrev[j] += fConnect[0] -> GetPreviousSolution(j);
            }
        }
    } else {
        for (int iside = 0; iside < compshape::NSides; iside++){
            
            int nshape = compshape::NShapeFunctions(iside,this->fMesh->GetDefaultOrder());
            if (nshape == 0) continue;
            int nstate = fMesh->NState();
            

            for (int ishape = 0; ishape < nshape; ishape++){
                double shapeFi = fIntegData.fPhi[count];
                for (int j = 0; j < nstate; j++ ){
                    fIntegData.fSol[j] += fConnect[iside] -> GetSolution(j) * shapeFi;
                    if (fIntegData.fSolPrev.size() != 0) {
                        fIntegData.fSolPrev[j] += fConnect[iside] -> GetPreviousSolution(j) * shapeFi;
                    }
                }
                count++;
            }
        }
    }

    // Store solution for nodes only in case the weakForm has no memory (used by topology optimization only) 
    if (fWeakForm && !fWeakForm->GetHasMemory()) {
        const size_t n_state = fMesh->NState();
        int nshape = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder());
        fIntegData.fSolNodes.resize(nshape * n_state);
        fIntegData.fSolNodes.setZero();
        for (int i = 0; i < nshape; i++){
            for (int j = 0; j < n_state; j++ ){
                fIntegData.fSolNodes[(n_state*i) + j] = fConnect[i]->GetSolution(j);
            }
        }
    }
}

template<class compshape>
void ElementT<compshape>::interpolateSolution(VecDouble &phi, VecDouble &u_) {
    u_.setZero();

    int count = 0;
    for (int iside = 0; iside < compshape::NSides; iside++){
        double shapeFi = phi[count];
        int nshape = compshape::NShapeFunctions(iside,this->fMesh->GetDefaultOrder());
        if (nshape == 0) continue;
        int nstate = fMesh->NState();

        for (int ishape = 0; ishape < nshape; ishape++){
            for (int j = 0; j < nstate; j++ ){
                u_[j] += fConnect[iside] -> GetSolution(j) * shapeFi;
            }
        }
        count += nshape;
    }
}

//------------------------------------------------------------------------------
//----------------------------INTERPOLATES VELOCITY-----------------------------
//------------------------------------------------------------------------------
template<class compshape>
void ElementT<compshape>::interpolateSolDerivatives(MatrixDouble &du_dx) {
    du_dx.setZero();    
    int DIM = compshape::Dimension;
    
    int count = 0;
    for (int iside = 0; iside < compshape::NSides; iside++){
        int nshape = compshape::NShapeFunctions(iside,this->fMesh->GetDefaultOrder());
        if (nshape == 0) continue;
        int nstate = fMesh->NState();

        for (int ishape = 0; ishape < nshape; ishape++){
            for (int j = 0; j < DIM; j++ ){
                for (int k = 0; k < nstate; k++ ){
                    du_dx(k,j) += fConnect[iside] -> GetSolution(k) *fIntegData.fDPhiX0(j,count);
                }
            }
        }
        count += nshape;
    }
    
}

template<class compshape>
void ElementT<compshape>::interpolateSolDerivatives(MatrixDouble &dphidx, MatrixDouble &du_dx) {
    du_dx.setZero();    
    int DIM = compshape::Dimension;
    int count = 0;
    for (int iside = 0; iside < compshape::NSides; iside++){
        int nshape = compshape::NShapeFunctions(iside,this->fMesh->GetDefaultOrder());
        if (nshape == 0) continue;
        int nstate = fMesh->NState();

        for (int ishape = 0; ishape < nshape; ishape++){
            for (int j = 0; j < DIM; j++ ){
                for (int k = 0; k < nstate; k++ ){
                    du_dx(k,j) += fConnect[iside] -> GetSolution(k) *dphidx(j,count);
                }
            }
        }
        count += nshape;
    }
}

template<class compshape>
void ElementT<compshape>::interpolateSolDerivatives() {
    fIntegData.fDSolDx.setZero();
    fIntegData.fDSolDAdim.setZero();
    bool flag = false;
    if (fIntegData.fDSolDxPrev.size() != 0){
        fIntegData.fDSolDxPrev.setZero(); 
        flag = true;
    }
    if (fConnect.size() == 0) return;

    if (fConnect.size() != compshape::NSides) {
        int nstate = fMesh->NState();
        for (int j = 0; j < nstate; j++ ){
            fIntegData.fDSolDx(j) = 0;
            if (flag) fIntegData.fDSolDxPrev(j) = 0;
            if (fIntegData.fNeedsDSolDAdim) fIntegData.fDSolDAdim(j) = 0;
        }
        return;
    } else {
        int count = 0;
        for (int iside = 0; iside < compshape::NSides; iside++){
            int nshape = compshape::NShapeFunctions(iside,this->fMesh->GetDefaultOrder());
            if (nshape == 0) continue;
            int nstate = fMesh->NState();

            for (int ishape = 0; ishape < nshape; ishape++){
                for (int j = 0; j < compshape::Dimension; j++ ){
                    for (int k = 0; k < nstate; k++ ){
                        fIntegData.fDSolDx(k,j) += fConnect[iside] -> GetSolution(k) *fIntegData.fDPhiX0(j,count);

                        if (flag) fIntegData.fDSolDxPrev(k,j) += fConnect[iside] -> GetPreviousSolution(k) * fIntegData.fDPhiX0(j,count);

                        if (fIntegData.fNeedsDSolDAdim) fIntegData.fDSolDAdim(k,j) += fConnect[iside] -> GetSolution(k) * fIntegData.fDPhi(j,count);
#ifdef DEBUG_BUILD
                        if (std::isnan(fIntegData.fDSolDx(k,j))){
                            PanicButton();
                        }
#endif
                    } 
                }
                count++;
            }
        }
    }

//     for (int i = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()); i--; ){
//         int nstate = fMesh->ConnectVec()[fConnect[i]]->GetNStateVariables();
//         for (int j = compshape::Dimension; j--; ){
//             for (int k = nstate; k--; ){
//                 fIntegData.fDSolDx(k,j) += fMesh->ConnectVec()[fConnect[i]] -> GetSolution(k) * fIntegData.fDPhiX0(j,i);
//                 if (flag) fIntegData.fDSolDxPrev(k,j) += fMesh->ConnectVec()[fConnect[i]] -> GetPreviousSolution(k) * fIntegData.fDPhiX0(j,i);

//                 if (fIntegData.fNeedsDSolDAdim) fIntegData.fDSolDAdim(k,j) += fMesh->ConnectVec()[fConnect[i]] -> GetSolution(k) * fIntegData.fDPhi(j,i);
// #ifdef DEBUG_BUILD
//                 if (std::isnan(fIntegData.fDSolDx(k,j))){
//                     PanicButton();
//                 }
// #endif
//             }
//         }
//     }
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//----------------MOUNTS EACH TYPE OF INCOMPRESSIBLE FLOW PROBEM----------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------- -----


//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<class compshape>
void ElementT<compshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

    if (!fWeakForm) return;

    // Already calculated before, re-use it.
    if(fIntegData.fStiffnessMatrix.has_value() && !fWeakForm->GetHasMemory()) {
        jacobianNRMatrix = *fIntegData.fStiffnessMatrix;
        rhsVector = fIntegData.fRHS;
        return;
    }

    int DIM = compshape::Dimension;
    fIntegData.fNeedsDSol = true;
    fIntegData.fDSolDx.resize(this->fWeakForm->NState(), DIM);
    fIntegData.fNeedsSol = true;
    fIntegData.fSol.resize(this->fWeakForm->NState());

    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (fWeakForm);
    auto *truss = dynamic_cast<PositionalTruss *> (fWeakForm);
    auto *frame2d = dynamic_cast<PositionalFrame2D *> (fWeakForm);

    if (frame2d){
        fIntegData.fDSolDAdim.resize(this->fWeakForm->NState(), DIM);
        fIntegData.fNeedsDSolDAdim = true;
    }
    
    int index = 0;
    double val = 0.;
    for(int it = 0; it < fIntRule.NPoints(); it++){
 
        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) fIntegData.fAdimCoord[k] = fIntRule.PointList(index,k);

        //Returns the quadrature integration weight
        fIntegData.fWeight = fIntRule.WeightList(index);

        //Computes the jacobian matrix
        fReference->ComputeJacobian(fIntegData);

        //Computes spatial derivatives
        ComputeSpatialDerivatives();
        
        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || truss || frame2d){
            fReference->ComputeCurrentJacobian(fIntegData,this);
            ComputeCurrentSpatialDerivatives();
        }

        if (fIntegData.fNeedsSol) interpolateSolution();
        if (fIntegData.fNeedsDSol) interpolateSolDerivatives();

        //Computes the element diffusion/viscosity matrix
        fWeakForm->ComputeStiffness(index, fIntegData, jacobianNRMatrix);
        
        //Computes the RHS vector
        fWeakForm->ComputeResidual(index, fIntegData, rhsVector); 

        index++;        
    };  

    // PrintMathematica(jacobianNRMatrix, "Stiffness");
    // PrintMathematica(rhsVector, "Rhs");
    // Set stiffness matrix to cache.
    // fIntegData.fStiffnessMatrix = jacobianNRMatrix;
    // fIntegData.fRHS = rhsVector;

    return;
};

template<class compshape>
void ElementT<compshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix){

    if (!fWeakForm) return;

    int DIM = compshape::Dimension;

    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (fWeakForm);
    auto *truss = dynamic_cast<PositionalTruss *> (fWeakForm);
    
    int index = 0;
    double val = 0.;
    for(int it = 0; it < fIntRule.NPoints(); it++){
 
        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) fIntegData.fAdimCoord[k] = fIntRule.PointList(index,k);

        //Returns the quadrature integration weight
        fIntegData.fWeight = fIntRule.WeightList(index);

        //Computes the jacobian matrix
        fReference->ComputeJacobian(fIntegData);

        //Computes spatial derivatives
        ComputeSpatialDerivatives();
        
        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || truss){
            fReference->ComputeCurrentJacobian(fIntegData,this);
            ComputeCurrentSpatialDerivatives();
        }
    
        //Computes the element diffusion/viscosity matrix
        fWeakForm->ComputeStiffness(index, fIntegData, jacobianNRMatrix);
        
        index++;        
    };  

    // std::cout << "Stiffness \n" << jacobianNRMatrix << std::endl;
    // std::cout << "Rhs \n" << rhsVector << std::endl;

    return;
};

template<class compshape>
void ElementT<compshape>::ComputeElContribution(VecDouble &rhsVector){

    if (!fWeakForm) return;

    int DIM = compshape::Dimension;
    fIntegData.fNeedsDSol = true;
    fIntegData.fDSolDx.resize(this->fWeakForm->NState(), DIM);
    fIntegData.fNeedsSol = true;
    fIntegData.fSol.resize(this->fWeakForm->NState());

    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (fWeakForm);
    auto *truss = dynamic_cast<PositionalTruss *> (fWeakForm);
    
    int index = 0;
    double val = 0.;
    for(int it = 0; it < fIntRule.NPoints(); it++){
 
        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) fIntegData.fAdimCoord[k] = fIntRule.PointList(index,k);

        //Returns the quadrature integration weight
        fIntegData.fWeight = fIntRule.WeightList(index);

        //Computes the jacobian matrix
        fReference->ComputeJacobian(fIntegData);

        //Computes spatial derivatives
        ComputeSpatialDerivatives();
        
        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || truss){
            fReference->ComputeCurrentJacobian(fIntegData,this);
            ComputeCurrentSpatialDerivatives();
        }     

        interpolateSolution();
        interpolateSolDerivatives();

        //Computes the RHS vector
        fWeakForm->ComputeResidual(index, fIntegData, rhsVector); 

        index++;        
    };  

    // std::cout << "Stiffness \n" << jacobianNRMatrix << std::endl;
    // std::cout << "Rhs \n" << rhsVector << std::endl;

    return;
};

//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<class compshape>
void ElementT<compshape>::ComputeElContribution(std::vector<MatrixDouble> &jacobianNRMatrix, std::vector<VecDouble> &rhsVector){

    if (!fWeakForm) return;

    int DIM = compshape::Dimension;
    fIntegData.fA0Inv.resize(DIM,DIM);
    fIntegData.fAdimCoord.resize(DIM);
    fIntegData.fNeedsDSol = true;
    fIntegData.fDSolDx.resize(this->fWeakForm->NState(), DIM);
    fIntegData.fNeedsSol = true;
    fIntegData.fSol.resize(this->fWeakForm->NState());

    int index = 0;
    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (fWeakForm);
    auto *truss = dynamic_cast<PositionalTruss *> (fWeakForm);

    for(int it = 0; it < fIntRule.NPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) fIntegData.fAdimCoord[k] = fIntRule.PointList(index,k);

        //Returns the quadrature integration weight
        fIntegData.fWeight = fIntRule.WeightList(index);

        //Computes the jacobian matrix
        fReference->ComputeJacobian(fIntegData);

        //Computes spatial derivatives
        ComputeSpatialDerivatives();
        
        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || truss){
            fReference->ComputeCurrentJacobian(fIntegData,this);
            ComputeCurrentSpatialDerivatives();
        }

        //Computes the element diffusion/viscosity matrix
        fWeakForm->ComputeStiffness(index, fIntegData, jacobianNRMatrix);

        if (fIntegData.fNeedsSol) interpolateSolution();
        if (fIntegData.fNeedsDSol) interpolateSolDerivatives();

        //Computes the RHS vector
        fWeakForm->ComputeResidual(index, fIntegData, rhsVector); 

        index++;        
    };  

    return;
};

//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<class compshape>
void ElementT<compshape>::ComputeElContribution(std::vector<MatrixDouble> &jacobianNRMatrix){

    if (!fWeakForm) return;

    int DIM = compshape::Dimension;
    fIntegData.fA0Inv.resize(DIM,DIM);
    fIntegData.fAdimCoord.resize(DIM);

    int index = 0;
    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (fWeakForm);
    auto *truss = dynamic_cast<PositionalTruss *> (fWeakForm);

    for(int it = 0; it < fIntRule.NPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) fIntegData.fAdimCoord[k] = fIntRule.PointList(index,k);

        //Returns the quadrature integration weight
        fIntegData.fWeight = fIntRule.WeightList(index);

        //Computes the jacobian matrix
        fReference->ComputeJacobian(fIntegData);

        //Computes spatial derivatives
        ComputeSpatialDerivatives();
        
        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || truss){
            fReference->ComputeCurrentJacobian(fIntegData,this);
            ComputeCurrentSpatialDerivatives();
        }

        //Computes the element diffusion/viscosity matrix
        fWeakForm->ComputeStiffness(index, fIntegData, jacobianNRMatrix); 

        index++;        
    };  
    return;
};

//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<class compshape>
void ElementT<compshape>::ComputeElContribution(std::vector<VecDouble> &rhsVector){

    if (!fWeakForm) return;

    int DIM = compshape::Dimension;
    fIntegData.fA0Inv.resize(DIM,DIM);
    fIntegData.fAdimCoord.resize(DIM);

    int index = 0;
    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (fWeakForm);
    auto *truss = dynamic_cast<PositionalTruss *> (fWeakForm);

    for(int it = 0; it < fIntRule.NPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) fIntegData.fAdimCoord[k] = fIntRule.PointList(index,k);

        //Returns the quadrature integration weight
        fIntegData.fWeight = fIntRule.WeightList(index);

        //Computes the jacobian matrix
        fReference->ComputeJacobian(fIntegData);

        //Computes spatial derivatives
        ComputeSpatialDerivatives();
        
        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || truss){
            fReference->ComputeCurrentJacobian(fIntegData,this);
            ComputeCurrentSpatialDerivatives();
        }

        if (fIntegData.fNeedsSol) interpolateSolution();
        if (fIntegData.fNeedsDSol) interpolateSolDerivatives();

        //Computes the RHS vector
        fWeakForm->ComputeResidual(index, fIntegData, rhsVector); 

        index++;        
    };  
    // // std::cout << "\nStiffness Element " << this->Index() << "\n" << jacobianNRMatrix[1];
    // // std::cout << "\nrhsVector Element " << this->Index() << "\n" << rhsVector[1];
   
    return;
};


//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<class compshape>
void ElementT<compshape>::ComputeError(VecDouble &errors){

#ifdef DEBUG_BUILD
    if (!fWeakForm || (fWeakForm->Dimension() != Mesh()->Dimension())) return;
    if (!fWeakForm->GetExactSolution()){
        std::cout << "Exact solution not set for material " << fWeakForm->Id() << std::endl;
        PanicButton();
    }
#endif

    int DIM = compshape::Dimension;
    fIntegData.fA0Inv.resize(DIM,DIM);
    fIntegData.fAdimCoord.resize(DIM);
    fIntegData.fSol.resize(fWeakForm->NState());
    fIntegData.fDSolDx.resize(fWeakForm->NState(),fWeakForm->Dimension());

    int index = 0;
    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (fWeakForm);
    auto *truss = dynamic_cast<PositionalTruss *> (fWeakForm);
    errors.setZero();

    for(int it = 0; it < fIntRule.NPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) fIntegData.fAdimCoord[k] = fIntRule.PointList(index,k);

        //Returns the quadrature integration weight
        fIntegData.fWeight = fIntRule.WeightList(index);

        //Computes the jacobian matrix
        fReference->ComputeJacobian(fIntegData);

        //Computes spatial derivatives
        ComputeSpatialDerivatives();
        
        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || truss){
            fReference->ComputeCurrentJacobian(fIntegData,this);
            ComputeCurrentSpatialDerivatives();
        }

        interpolateSolution();
        interpolateSolDerivatives();

        //Computes the RHS vector
        fWeakForm->ComputeError(fIntegData, errors); 

        index++;        
    };  
    return;
};

template<class compshape >
Element * ElementT<compshape>::Clone() const {
    return new ElementT(*this);
}

template<class compshape >
void ElementT<compshape>::Integrate(std::vector<std::string> &varNames, std::map<std::string,VecDouble> &result){
    int DIM = compshape::Dimension;
    fIntegData.fA0Inv.resize(DIM,DIM);
    fIntegData.fAdimCoord.resize(DIM);
    fIntegData.fSol.resize(fWeakForm->NState());
    fIntegData.fDSolDx.resize(fWeakForm->NState(),fWeakForm->Dimension());

    
    for (int ivar = 0; ivar < varNames.size(); ivar++){
        int varindex = fWeakForm->VariableIndex(varNames[ivar]);
        int nvar = fWeakForm->NSolutionVariables(varindex);
        VecDouble Sol(nvar);
        Sol.setZero();
        if (result[varNames[ivar]].size() == 0) result[varNames[ivar]] = Sol;
        int index = 0;
        for(int it = 0; it < fIntRule.NPoints(); it++){

            //Defines the integration points adimentional coordinates
            for (int k = 0; k < compshape::Dimension; k++) fIntegData.fAdimCoord[k] = fIntRule.PointList(index,k);

            //Returns the quadrature integration weight
            fIntegData.fWeight = fIntRule.WeightList(index);

            //Computes the jacobian matrix
            fReference->ComputeJacobian(fIntegData);

            //Computes spatial derivatives
            ComputeSpatialDerivatives();
            
            // Computes current spatial derivatives (only for position-based weak forms)
            // if (pos2d || truss){
            //     fReference->ComputeCurrentJacobian(fIntegData,this);
            //     ComputeCurrentSpatialDerivatives();
            // }

            interpolateSolution();
            interpolateSolDerivatives();

            double WJ = fIntegData.fWeight * fIntegData.fJacA0;
            Sol.setZero();
            fIntegData.fIndex = index;
            fWeakForm->Solution(fIntegData,varindex,Sol);
            // Solution(varindex,Sol);
            result[varNames[ivar]] += Sol * WJ;
            index++;
        }
    }
}

#include "ShapeHexahedron.h"
#include "ShapeOneDLin.h"
#include "ShapeOneDQua.h"
#include "ShapeOneDCub.h"
#include "ShapeQuadrilateralLin.h"
#include "ShapeQuadrilateralQua.h"
#include "ShapePoint.h"
#include "ShapeTetrahedronLin.h"
#include "ShapeTetrahedronQua.h"
#include "ShapeTetrahedronCub.h"
#include "ShapeTriangleLin.h"
#include "ShapeTriangleQua.h"
#include "ShapeTriangleCub.h"
#include "HierarquicalOneD.h"
#include "HierarquicalQuad.h"
#include "HierarquicalTriangle.h"

template class ElementT<ShapePoint>;
template class ElementT<ShapeOneDLin>;
template class ElementT<ShapeOneDQua>;
template class ElementT<ShapeOneDCub>;
template class ElementT<HierarquicalOneD>;
template class ElementT<ShapeTriangleLin>;
template class ElementT<ShapeTriangleQua>;
template class ElementT<ShapeTriangleCub>;
template class ElementT<HierarquicalTriangle>;
template class ElementT<ShapeQuadrilateralLin>;
template class ElementT<ShapeQuadrilateralQua>;
template class ElementT<HierarquicalQuad>;
template class ElementT<ShapeTetrahedronLin>;
template class ElementT<ShapeTetrahedronQua>;
template class ElementT<ShapeTetrahedronCub>;
template class ElementT<ShapeHexahedron>;