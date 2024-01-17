#include "ElementT.h"
#include "PositionalTruss.h"
#include "ElasticityPositional2D.h"

// #include "Boundary.h"
// #include "CompMesh.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

template<class tshape>
void ElementT<tshape>::setIntegPointWeightFunction() {
    //It produces wrong results for constant weight functions
    VecDouble xsi(tshape::Dimension);    

    for(int i = 0; i < fIntRule.NPoints(); i++) {
        fIntegData.fPrevWeightFunction[i] = fIntegData.fWeightFunction[i];
        fIntegData.fWeightFunction[i] = 0.;
    }

    int index=0;

    for(int it = 0; it < fIntRule.NPoints(); it++){
        
       xsi[0] = fIntRule.PointList(index,0);
       xsi[1] = fIntRule.PointList(index,1);

       for (int j=0; j<tshape::NElNodes; j++){
           fIntegData.fWeightFunction[index] += fIntegData.fPhi[j] * fMesh->NodeVec()[fConnect[j]] -> getWeightFunction();
       };
       // fIntegData.fWeightFunction(index) = 1.;
       index++;
    }; 

    // index = 0;
     return;
};

//------------------------------------------------------------------------------
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------


template<class tshape>
void ElementT<tshape>::ComputeIntPointDistFunction(VecDouble &nodalval) {
    
    int DIM = tshape::Dimension;
    fIntegData.fAdimCoord.resize(DIM);
    fIntegData.fPhi.resize(tshape::NElNodes);
    
    // for(int i = 0; i < nQuad.getNumberOfIntegrationPoints(); i++) {
    //     intPointWeightFunctionPrev[i] = intPointWeightFunction[i];
    //     intPointWeightFunction[i] = 0.;
    // }

    int index=0;
    fIntegData.fDistFunction.setZero();
    for(int it = 0; it < fIntRule.NPoints(); it++){
        
        for (int i=0; i<DIM; i++) fIntegData.fAdimCoord[i] = fIntRule.PointList(index,i);
           
        tshape::Shape(fIntegData.fAdimCoord,fIntegData.fPhi);

        for (int j=0; j<tshape::NElNodes; j++){
                fIntegData.fDistFunction[index] +=  fIntegData.fPhi[j] * nodalval[j];
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

    fIntPointCoordinates.resize(fIntRule.NPoints(),2);
    fIntPointCoordinates.setZero();
    
    VecDouble xsi(DIM);
    VecDouble phi_(tshape::NElNodes);
    fIntPointCoordinates.resize(fIntRule.NPoints(),DIM);

    for (int i = 0; i < fIntRule.NPoints(); i++){
        double x[DIM] = {};

        for (int k = DIM; k--; ) xsi[k] = fIntRule.PointList(i,k);

        tshape::Shape(xsi,phi_);

        for (int k = DIM; k--; ) fIntPointCoordinates(i,k) = 0.;

        for (int j = 0; j < tshape::NElNodes; j++)
            for (int k = DIM; k--; )
                fIntPointCoordinates(i,k) += fMesh->NodeVec()[fConnect[j]] -> getCoordinateValue(k) * phi_[j];
        
    };

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
void ElementT<tshape>::ComputeJacobian() {

    int DIM = tshape::Dimension;
    fIntegData.fA0Inv.resize(DIM,DIM);
    fIntegData.fA0Inv.setZero();
    fIntegData.fAxes0.resize(3,DIM);
    fIntegData.fAxes0.setZero();
    fIntegData.fA0.resize(DIM,DIM);
    fIntegData.fA0.setZero();
    fIntegData.fX.resize(3);
    fIntegData.fX.setZero();

    fIntegData.fPhi.resize(tshape::NElNodes);
    fIntegData.fPhi.setZero();

    fIntegData.fDPhi.resize(tshape::NElNodes,tshape::Dimension);
    fIntegData.fDPhi.setZero();

    tshape::Shape(fIntegData.fAdimCoord,fIntegData.fPhi);
    tshape::ShapeGradient(fIntegData.fAdimCoord,fIntegData.fDPhi);   
   
    fIntegData.fA0.setZero();
    VecDouble xna(3);
    for (int i = tshape::NElNodes; i--; ){
        for (int j = 3; j--; ){
            // Approximate the integration space
            fIntegData.fX[j] += fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j) * fIntegData.fPhi(i);
            xna[j] = fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j);
            
            for (int k = DIM; k--; ){
                fIntegData.fAxes0(j,k) += xna[j] * fIntegData.fDPhi(i,k);
            };
        };
    };

    int ncols = DIM;
    int dim = DIM;

    switch (DIM) {
        case 0:
            fIntegData.fJacA0 = 1.;
            break;
        case 1:
        {
            VecDouble v_1(3);
            v_1.setZero();

            for (int i = 0; i < 3; i++) {
                v_1[i] = fIntegData.fAxes0(i, 0);
            }

            double norm_v_1 = 0.;
            for (int i = 0; i < 3; i++) {
                norm_v_1 += v_1[i] * v_1[i];
            }

            norm_v_1 = sqrt(norm_v_1);
            fIntegData.fA0(0, 0) = norm_v_1;
            fIntegData.fJacA0 = norm_v_1;
            fIntegData.fA0Inv(0, 0) = 1.0 / fIntegData.fJacA0;

            fIntegData.fJacA0 = fabs(fIntegData.fJacA0);

        }
            break;
        case 2:
        {
            //     //Computing the jacobian determinant and Inverse
            fIntegData.fA0(0,0) = fIntegData.fAxes0(0,0);
            fIntegData.fA0(0,1) = fIntegData.fAxes0(0,1);
            fIntegData.fA0(1,0) = fIntegData.fAxes0(1,0);
            fIntegData.fA0(1,1) = fIntegData.fAxes0(1,1);
            fIntegData.fJacA0 = fIntegData.fA0(0,0) * fIntegData.fA0(1,1) - fIntegData.fA0(0,1) * fIntegData.fA0(1,0);

            fIntegData.fA0Inv(0,0) = fIntegData.fA0(1,1) / fIntegData.fJacA0;
            fIntegData.fA0Inv(1,1) = fIntegData.fA0(0,0) / fIntegData.fJacA0;
            fIntegData.fA0Inv(0,1) = -fIntegData.fA0(0,1) / fIntegData.fJacA0;
            fIntegData.fA0Inv(1,0) = -fIntegData.fA0(1,0) / fIntegData.fJacA0;
            fIntegData.fJacA0 = fabs(fIntegData.fJacA0);
        }
            break;
        case 3:
        {
            // axes.resize(dim, 3);

            for (int i = 0; i < 3; i++) {
                fIntegData.fA0(i, 0) = fIntegData.fAxes0(i, 0);
                fIntegData.fA0(i, 1) = fIntegData.fAxes0(i, 1);
                fIntegData.fA0(i, 2) = fIntegData.fAxes0(i, 2);
            }

            fIntegData.fJacA0 -= fIntegData.fA0(0, 2) * fIntegData.fA0(1, 1) * fIntegData.fA0(2, 0); //- a02 a11 a20
            fIntegData.fJacA0 += fIntegData.fA0(0, 1) * fIntegData.fA0(1, 2) * fIntegData.fA0(2, 0); //+ a01 a12 a20
            fIntegData.fJacA0 += fIntegData.fA0(0, 2) * fIntegData.fA0(1, 0) * fIntegData.fA0(2, 1); //+ a02 a10 a21
            fIntegData.fJacA0 -= fIntegData.fA0(0, 0) * fIntegData.fA0(1, 2) * fIntegData.fA0(2, 1); //- a00 a12 a21
            fIntegData.fJacA0 -= fIntegData.fA0(0, 1) * fIntegData.fA0(1, 0) * fIntegData.fA0(2, 2); //- a01 a10 a22
            fIntegData.fJacA0 += fIntegData.fA0(0, 0) * fIntegData.fA0(1, 1) * fIntegData.fA0(2, 2); //+ a00 a11 a22

            fIntegData.fA0Inv(0, 0) = (-fIntegData.fA0(1, 2) * fIntegData.fA0(2, 1) + fIntegData.fA0(1, 1) * fIntegData.fA0(2, 2)) / fIntegData.fJacA0; //-a12 a21 + a11 a22
            fIntegData.fA0Inv(0, 1) = (fIntegData.fA0(0, 2) * fIntegData.fA0(2, 1) - fIntegData.fA0(0, 1) * fIntegData.fA0(2, 2)) / fIntegData.fJacA0; //a02 a21 - a01 a22
            fIntegData.fA0Inv(0, 2) = (-fIntegData.fA0(0, 2) * fIntegData.fA0(1, 1) + fIntegData.fA0(0, 1) * fIntegData.fA0(1, 2)) / fIntegData.fJacA0; //-a02 a11 + a01 a12
            fIntegData.fA0Inv(1, 0) = (fIntegData.fA0(1, 2) * fIntegData.fA0(2, 0) - fIntegData.fA0(1, 0) * fIntegData.fA0(2, 2)) / fIntegData.fJacA0; //a12 a20 - a10 a22
            fIntegData.fA0Inv(1, 1) = (-fIntegData.fA0(0, 2) * fIntegData.fA0(2, 0) + fIntegData.fA0(0, 0) * fIntegData.fA0(2, 2)) / fIntegData.fJacA0; //-a02 a20 + a00 a22
            fIntegData.fA0Inv(1, 2) = (fIntegData.fA0(0, 2) * fIntegData.fA0(1, 0) - fIntegData.fA0(0, 0) * fIntegData.fA0(1, 2)) / fIntegData.fJacA0; //a02 a10 - a00 a12
            fIntegData.fA0Inv(2, 0) = (-fIntegData.fA0(1, 1) * fIntegData.fA0(2, 0) + fIntegData.fA0(1, 0) * fIntegData.fA0(2, 1)) / fIntegData.fJacA0; //-a11 a20 + a10 a21
            fIntegData.fA0Inv(2, 1) = (fIntegData.fA0(0, 1) * fIntegData.fA0(2, 0) - fIntegData.fA0(0, 0) * fIntegData.fA0(2, 1)) / fIntegData.fJacA0; //a01 a20 - a00 a21
            fIntegData.fA0Inv(2, 2) = (-fIntegData.fA0(0, 1) * fIntegData.fA0(1, 0) + fIntegData.fA0(0, 0) * fIntegData.fA0(1, 1)) / fIntegData.fJacA0; //-a01 a10 + a00 a11

            fIntegData.fJacA0 = fabs(fIntegData.fJacA0);

            // axes.setZero();
            // axes(0, 0) = 1.0;
            // axes(1, 1) = 1.0;
            // axes(2, 2) = 1.0;
        }
            break;
    }

    return;
};

template<class tshape>
void ElementT<tshape>::ComputeCurrentJacobian() {

    int DIM = tshape::Dimension;
    fIntegData.fA1.resize(DIM,DIM);
    fIntegData.fA1.setZero();
    fIntegData.fX1.resize(3);
    fIntegData.fX1.setZero();
    fIntegData.fAxes1.resize(3,DIM);
    fIntegData.fAxes1.setZero();

    VecDouble yna(3);
    
    fIntegData.fA1.setZero();
    for (int i = tshape::NElNodes; i--; ){
        for (int j = fWeakForm->NState(); j--; ){
            // Approximate the integration space
            fIntegData.fX1[j] += (fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j) + fMesh->NodeVec()[fConnect[i]] -> GetSolution(j)) * fIntegData.fPhi(i);
            yna[j] = fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j) + fMesh->NodeVec()[fConnect[i]] -> GetSolution(j);

            for (int k = DIM; k--; ){
                fIntegData.fAxes1(j,k) += yna[j] * fIntegData.fDPhi(i,k);
            };
        };
    };

    switch (tshape::Dimension)
    {
    case 1:
    {
        VecDouble v_1(3);
        v_1.setZero();

        for (int i = 0; i < 3; i++) {
            v_1[i] = fIntegData.fAxes1(i, 0);
        }

        double norm_v_1 = 0.;
        for (int i = 0; i < 3; i++) {
            norm_v_1 += v_1[i] * v_1[i];
        }
        norm_v_1 = sqrt(norm_v_1);
        fIntegData.fA1(0, 0) = norm_v_1;
        fIntegData.fJacA1 = norm_v_1;
        fIntegData.fJacA1 = fabs(fIntegData.fJacA1);
    }
        break;
    case 2:
    {
        //Computing the jacobian determinant and Inverse
        fIntegData.fA1(0,0) = fIntegData.fAxes1(0,0);
        fIntegData.fA1(0,1) = fIntegData.fAxes1(0,1);
        fIntegData.fA1(1,0) = fIntegData.fAxes1(1,0);
        fIntegData.fA1(1,1) = fIntegData.fAxes1(1,1);
        fIntegData.fJacA1 = fIntegData.fA1(0,0) * fIntegData.fA1(1,1) - fIntegData.fA1(0,1) * fIntegData.fA1(1,0);
        fIntegData.fJacA1 = fabs(fIntegData.fJacA1);
    }
        break;
    case 3:
    {
        // axes.resize(dim, 3);

        for (int i = 0; i < 3; i++) {
            fIntegData.fA1(i, 0) = fIntegData.fAxes1(i, 0);
            fIntegData.fA1(i, 1) = fIntegData.fAxes1(i, 1);
            fIntegData.fA1(i, 2) = fIntegData.fAxes1(i, 2);
        }

        fIntegData.fJacA1 -= fIntegData.fA1(0, 2) * fIntegData.fA1(1, 1) * fIntegData.fA1(2, 0); //- a02 a11 a20
        fIntegData.fJacA1 += fIntegData.fA1(0, 1) * fIntegData.fA1(1, 2) * fIntegData.fA1(2, 0); //+ a01 a12 a20
        fIntegData.fJacA1 += fIntegData.fA1(0, 2) * fIntegData.fA1(1, 0) * fIntegData.fA1(2, 1); //+ a02 a10 a21
        fIntegData.fJacA1 -= fIntegData.fA1(0, 0) * fIntegData.fA1(1, 2) * fIntegData.fA1(2, 1); //- a00 a12 a21
        fIntegData.fJacA1 -= fIntegData.fA1(0, 1) * fIntegData.fA1(1, 0) * fIntegData.fA1(2, 2); //- a01 a10 a22
        fIntegData.fJacA1 += fIntegData.fA1(0, 0) * fIntegData.fA1(1, 1) * fIntegData.fA1(2, 2); //+ a00 a11 a22

        fIntegData.fJacA1 = fabs(fIntegData.fJacA1);

        // axes.setZero();
        // axes(0, 0) = 1.0;
        // axes(1, 1) = 1.0;
        // axes(2, 2) = 1.0;
    }
    
    default:
        std::cout << "Please implement me";
        PanicButton();
        break;
    }
    return;
};

//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::ComputeSpatialDerivatives() {
    
    
        
    tshape::ShapeGradient(fIntegData.fAdimCoord,fIntegData.fDPhi);
    // shapeQuad.ShapeHessian(xsi,ddphi);
    
    fIntegData.fDPhiX0.setZero();

    //Shape functions spatial first derivatives
    fIntegData.fDPhiX0 = fIntegData.fDPhi * fIntegData.fA0Inv;

    return;
};

template<class tshape>
void ElementT<tshape>::ComputeCurrentSpatialDerivatives() {
    
    fIntegData.fDPhi.resize(tshape::NElNodes,tshape::Dimension);
    fIntegData.fDPhi.setZero();
    
    tshape::ShapeGradient(fIntegData.fAdimCoord,fIntegData.fDPhi);
    // shapeQuad.ShapeHessian(xsi,ddphi);
    
    fIntegData.fDPhiX1.setZero();

    //Shape functions spatial first derivatives
    fIntegData.fDPhiX1 = fIntegData.fDPhi * fIntegData.fA1.inverse().transpose();

    return;
};
//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::getHighOrderSpatialDerivatives(VecDouble &xsi, MatrixDouble &ainv_, MatrixDouble &dphi_dx, MatrixDouble &dDphi_dx) {
    
    int DIM = tshape::Dimension;
    dDphi_dx.setZero();
    std::vector<MatrixDouble> ddphi(tshape::NElNodes,MatrixDouble(DIM,DIM));
    for (int i = 0; i < tshape::NElNodes; i++)
    {
        ddphi[i].setZero();
    }
    
    tshape::ShapeHessian(xsi,ddphi);
    
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
    VecDouble xna_(DIM);

    for (int i = tshape::NElNodes; i--; ){
        xna_.setZero();
        for (int j = DIM; j--; ){
            // Approximate the integration space
            xna_[j] = fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j);
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
        for (int i = tshape::NElNodes; i--; ){
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
    
    for (int i = tshape::NElNodes; i--; ){
        double shapeFi = fIntegData.fPhi[i];
        for (int j = DIM; j--; ){
            umesh_[j] += fMesh->NodeVec()[fConnect[i]] -> getMeshVelocity(j) * shapeFi;
        }
    }

    return;
}

template<class tshape>
void ElementT<tshape>::interpolateSolution(int &index, VecDouble &u_) {
    u_.setZero();
    for (int i = tshape::NElNodes; i--; ){
        double shapeFi = fIntegData.fPhi[i];
        int nstate = fMesh->NodeVec()[fConnect[i]]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            u_[j] += fMesh->NodeVec()[fConnect[i]] -> GetSolution(j) * shapeFi;
        }
    }
}
template<class tshape>
double ElementT<tshape>::InterpolateVariable(VecDouble &nValues, int point) {
    double val = 0.;
    fIntegData.fPhi.resize(tshape::NElNodes);
    fIntegData.fPhi.setZero();
    fIntegData.fAdimCoord.resize(tshape::Dimension);
    for (int i=0; i<tshape::Dimension; i++) fIntegData.fAdimCoord[i] = fIntRule.PointList(point,i);
    tshape::Shape(fIntegData.fAdimCoord,fIntegData.fPhi);
    for (int i = tshape::NElNodes; i--; ){
        double shapeFi = fIntegData.fPhi[i];
        val += nValues[i] * shapeFi;
    }
    return val;
}

template<class tshape>
void ElementT<tshape>::interpolateSolDTimeDerivatives() {
    fIntegData.fDSolDt.setZero();
    fIntegData.fDSolDDt.setZero();
    for (int i = tshape::NElNodes; i--; ){
        double shapeFi = fIntegData.fPhi[i];
        int nstate = fMesh->NodeVec()[fConnect[i]]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            fIntegData.fDSolDt[j] += fMesh->NodeVec()[fConnect[i]] -> GetDSolutionDTime(j) * shapeFi;
            fIntegData.fDSolDDt[j] += fMesh->NodeVec()[fConnect[i]] -> GetDSolutionDDTime(j) * shapeFi;
        }
    }
}

template<class tshape>
void ElementT<tshape>::interpolateSolDTimeDerivatives(VecDouble &du_dt, VecDouble &du_ddt) {
    fIntegData.fDSolDt.setZero();
    fIntegData.fDSolDDt.setZero();
    for (int i = tshape::NElNodes; i--; ){
        double shapeFi = fIntegData.fPhi[i];
        int nstate = fMesh->NodeVec()[fConnect[i]]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            fIntegData.fDSolDt[j] += fMesh->NodeVec()[fConnect[i]] -> GetDSolutionDTime(j) * shapeFi;
            fIntegData.fDSolDDt[j] += fMesh->NodeVec()[fConnect[i]] -> GetDSolutionDDTime(j) * shapeFi;
        }
    }
    du_dt = fIntegData.fDSolDt;
    du_ddt = fIntegData.fDSolDDt;
}


template<class tshape>
void ElementT<tshape>::interpolateSolution() {
    fIntegData.fSol.setZero();
    for (int i = tshape::NElNodes; i--; ){
        double shapeFi = fIntegData.fPhi[i];
        int nstate = fMesh->NodeVec()[fConnect[i]]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            fIntegData.fSol[j] += fMesh->NodeVec()[fConnect[i]] -> GetSolution(j) * shapeFi;
        }
    }
}

template<class tshape>
void ElementT<tshape>::interpolateSolution(VecDouble &phi, VecDouble &u_) {
    u_.setZero();
    for (int i = tshape::NElNodes; i--; ){
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
    for (int i = tshape::NElNodes; i--; ){
        int nstate = fMesh->NodeVec()[fConnect[i]]->GetNStateVariables();
        for (int j = DIM; j--; ){
            for (int k = nstate; k--; ){
                du_dx(k,j) += fMesh->NodeVec()[fConnect[i]] -> GetSolution(k) * fIntegData.fDPhiX0(i,j);
            }
        }
    }
}

template<class tshape>
void ElementT<tshape>::interpolateSolDerivatives(MatrixDouble &dphidx, MatrixDouble &du_dx) {
    du_dx.setZero();    
    int DIM = tshape::Dimension;
    for (int i = tshape::NElNodes; i--; ){
        int nstate = fMesh->NodeVec()[fConnect[i]]->GetNStateVariables();
        for (int j = DIM; j--; ){
            for (int k = nstate; k--; ){
                du_dx(k,j) += fMesh->NodeVec()[fConnect[i]] -> GetSolution(k) * dphidx(i,j);
            }
        }
    }
}

template<class tshape>
void ElementT<tshape>::interpolateSolDerivatives() {
    fIntegData.fDSolDx.setZero();    
    int DIM = tshape::Dimension;
    for (int i = tshape::NElNodes; i--; ){
        int nstate = fMesh->NodeVec()[fConnect[i]]->GetNStateVariables();
        for (int j = DIM; j--; ){
            for (int k = nstate; k--; ){
                fIntegData.fDSolDx(k,j) += fMesh->NodeVec()[fConnect[i]] -> GetSolution(k) * fIntegData.fDPhiX0(i,j);
            }
        }
    }
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
template<class tshape>
void ElementT<tshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

    if (!fWeakForm) return;

    int DIM = tshape::Dimension;

    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (fWeakForm);
    auto *truss = dynamic_cast<PositionalTruss *> (fWeakForm);
    
    int index = 0;
    fIntegData.fAdimCoord.resize(DIM);
    for(int it = 0; it < fIntRule.NPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) fIntegData.fAdimCoord[k] = fIntRule.PointList(index,k);

        //Returns the quadrature integration weight
        fIntegData.fWeight = fIntRule.WeightList(index);

        //Computes the jacobian matrix
        ComputeJacobian();

        //Computes spatial derivatives
        ComputeSpatialDerivatives();
        
        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || truss){
            ComputeCurrentJacobian();
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

    // std::cout << "Stiffness \n" << jacobianNRMatrix << std::endl;
    // std::cout << "Rhs \n" << rhsVector << std::endl;

    return;
};


//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::ComputeElContribution(std::vector<MatrixDouble> &jacobianNRMatrix, std::vector<VecDouble> &rhsVector){

    if (!fWeakForm) return;

    int DIM = tshape::Dimension;
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
        ComputeJacobian();

        //Computes spatial derivatives
        ComputeSpatialDerivatives();
        
        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || truss){
            ComputeCurrentJacobian();
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
    // // std::cout << "\nStiffness Element " << this->Index() << "\n" << jacobianNRMatrix[1];
    // // std::cout << "\nrhsVector Element " << this->Index() << "\n" << rhsVector[1];
    // //Apply boundary conditions
    // ApplyBC(jacobianNRMatrix, rhsVector);

    return;
};


//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementT<tshape>::ComputeError(VecDouble &errors){

    if (!fWeakForm || (fWeakForm->Dimension() != Mesh()->Dimension())) return;
    if (!fWeakForm->GetExactSolution()){
        std::cout << "Exact solution not set for material " << fWeakForm->Id() << std::endl;
        PanicButton();
    }

    int DIM = tshape::Dimension;
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
        ComputeJacobian();

        //Computes spatial derivatives
        ComputeSpatialDerivatives();
        
        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || truss){
            ComputeCurrentJacobian();
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

template<class tshape >
Element * ElementT<tshape>::Clone() const {
    return new ElementT(*this);
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

template class ElementT<ShapePoint>;
template class ElementT<ShapeOneDLin>;
template class ElementT<ShapeOneDQua>;
template class ElementT<ShapeOneDCub>;
template class ElementT<ShapeTriangleLin>;
template class ElementT<ShapeTriangleQua>;
template class ElementT<ShapeTriangleCub>;
template class ElementT<ShapeQuadrilateralLin>;
template class ElementT<ShapeQuadrilateralQua>;
template class ElementT<ShapeTetrahedronLin>;
template class ElementT<ShapeTetrahedronQua>;
template class ElementT<ShapeTetrahedronCub>;
template class ElementT<ShapeHexahedron>;