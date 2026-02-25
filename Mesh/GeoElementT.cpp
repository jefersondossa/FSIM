#include "GeoElementT.h"

template<class geoshape>
GeoElementT<geoshape>::GeoElementT() : GeoElement(){

};

template<class geoshape>
GeoElementT<geoshape>::GeoElementT(int64_t index, VecInt &geonodes, GeoMesh* mesh, int &material) : GeoElement(){
    fMesh = mesh;
    fGeoNodes.resize(geoshape::NCornerNodes);
    fIndex = index;
    for (int i = geoshape::NCornerNodes; i--; ) fGeoNodes[i] = geonodes[i];
    fMaterial = material;
    fNeighborElements.clear();
};



//------------------------------------------------------------------------------
//-------------------------SPATIAL TRANSFORM - JACOBIAN-------------------------
//------------------------------------------------------------------------------
template<class geohape>
void GeoElementT<geoshape>::ComputeJacobian() {

    int DIM = geoshape::Dimension;
    fIntegData.fA0Inv.setZero();
    fIntegData.fAxes0.setZero();
    fIntegData.fA0.setZero();
    fIntegData.fX.setZero();
    VecDouble phigeo(geoshape::NShape);
    MatrixDouble dphigeo(geoshape::Dimension,geoshape::NShape);
    VecInt orders(geoshape::NSides);
    orders.fill(1);

    geoshape::Shape(fIntegData.fAdimCoord,phigeo,dphigeo,orders);

    fIntegData.fA0.setZero();
    VecDouble xna(3);
    for (int i = geoshape::NShape; i--; ){
        for (int j = 3; j--; ){
            // Approximate the integration space
            fIntegData.fX[j] += fMesh->Reference()->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) * phigeo(i);
            // std::cout << "Coord " << i << " " << j << " = " << fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j) << std::endl;
            xna[j] = fMesh->Reference()->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j);
            
            for (int k = DIM; k--; ){
                fIntegData.fAxes0(j,k) += xna[j] * dphigeo(k,i);
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
            auto *BC = dynamic_cast<L2Projection *> (fWeakForm);

            if (BC){
                VecDouble v_1(3), v_2(3);
                VecDouble v_1_til(3), v_2_til(3);
                v_1.setZero();
                v_2.setZero();
                v_1_til.setZero();
                v_2_til.setZero();

                for (int i = 0; i < 3; i++) {
                    v_1[i] = fIntegData.fAxes0(i, 0);
                    v_2[i] = fIntegData.fAxes0(i, 1);
                    // v_1[i] = gradx(i, 0);
                    // v_2[i] = gradx(i, 1);
                }

                double norm_v_1_til = 0.0;
                double norm_v_2_til = 0.0;
                double v_1_dot_v_2 = 0.0;

                for (int i = 0; i < 3; i++) {
                    norm_v_1_til += v_1[i] * v_1[i];
                    v_1_dot_v_2 += v_1[i] * v_2[i];
                }
                norm_v_1_til = sqrt(norm_v_1_til);

                for (int i = 0; i < 3; i++) {
                    v_1_til[i] = v_1[i] / norm_v_1_til;
                    v_2_til[i] = v_2[i] - v_1_dot_v_2 * v_1_til[i] / norm_v_1_til;
                    norm_v_2_til += v_2_til[i] * v_2_til[i];
                }
                norm_v_2_til = sqrt(norm_v_2_til);


                fIntegData.fA0(0, 0) = norm_v_1_til;
                fIntegData.fA0(0, 1) = v_1_dot_v_2 / norm_v_1_til;
                fIntegData.fA0(1, 1) = norm_v_2_til;

                fIntegData.fJacA0 = fIntegData.fA0(0, 0) * fIntegData.fA0(1, 1) - fIntegData.fA0(1, 0) * fIntegData.fA0(0, 1);

                fIntegData.fA0Inv(0, 0) = +fIntegData.fA0(1, 1) / fIntegData.fJacA0;
                fIntegData.fA0Inv(1, 1) = +fIntegData.fA0(0, 0) / fIntegData.fJacA0;
                fIntegData.fA0Inv(0, 1) = -fIntegData.fA0(0, 1) / fIntegData.fJacA0;
                fIntegData.fA0Inv(1, 0) = -fIntegData.fA0(1, 0) / fIntegData.fJacA0;

                fIntegData.fJacA0 = fabs(fIntegData.fJacA0);

                for (int i = 0; i < 3; i++) {
                    v_2_til[i] /= norm_v_2_til;
                    fIntegData.fAxes0(i,0) = v_1_til[i];
                    fIntegData.fAxes0(i,1) = v_2_til[i];
                }

                // std::cout << "fIntegData.fJacA0 = " << fIntegData.fJacA0 << std::endl;
                // std::cout << "fIntegData.fA0 = \n" << fIntegData.fA0 << std::endl;
                // std::cout << "fIntegData.fA0Inv = \n" << fIntegData.fA0Inv << std::endl;
            } else {
                //Computing the jacobian determinant and Inverse
                fIntegData.fA0(0,0) = fIntegData.fAxes0(0,0);
                fIntegData.fA0(0,1) = fIntegData.fAxes0(0,1);
                fIntegData.fA0(1,0) = fIntegData.fAxes0(1,0);
                fIntegData.fA0(1,1) = fIntegData.fAxes0(1,1);
                fIntegData.fJacA0 = fIntegData.fA0(0,0) * fIntegData.fA0(1,1) - fIntegData.fA0(0,1) * fIntegData.fA0(1,0);

                fIntegData.fA0Inv(0,0) = fIntegData.fA0(1,1) / fIntegData.fJacA0;
                fIntegData.fA0Inv(1,1) = fIntegData.fA0(0,0) / fIntegData.fJacA0;
                fIntegData.fA0Inv(0,1) = -fIntegData.fA0(0,1) / fIntegData.fJacA0;
                fIntegData.fA0Inv(1,0) = -fIntegData.fA0(1,0) / fIntegData.fJacA0;

                // std::cout << "fIntegData.fJacA0 = " << fIntegData.fJacA0 << std::endl;
                // std::cout << "fIntegData.fA0 = \n" << fIntegData.fA0 << std::endl;
                // std::cout << "fIntegData.fA0Inv = \n" << fIntegData.fA0Inv << std::endl;

                fIntegData.fJacA0 = fabs(fIntegData.fJacA0);
            }
            
#ifdef DEBUG_BUILD
            if (fIntegData.fJacA0 < 1.e-6) PanicButton();
#endif
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
            fIntegData.fJacA0 = 0.;
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


template<class geohape>
void GeoElementT<geoshape>::ComputeJacobianSearch() {

    int DIM = geoshape::Dimension;
    fIntegData.fA0Inv.setZero();
    fIntegData.fAxes0.setZero();
    fIntegData.fA0.setZero();
    fIntegData.fX.setZero();
    fIntegData.fPhi.setZero();
    fIntegData.fDPhi.setZero();
    VecDouble phigeo(geoshape::NShape);
    MatrixDouble dphigeo(geoshape::Dimension,geoshape::NShape);
    VecInt orders(geoshape::NSides);
    orders.fill(1);

    geoshape::Shape(fIntegData.fAdimCoord,phigeo,dphigeo,orders);

    fIntegData.fA0.setZero();
    VecDouble xna(3);
    for (int i = geoshape::NShape; i--; ){
        for (int j = 3; j--; ){
            // Approximate the integration space
            fIntegData.fX[j] += fMesh->Reference()->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) * phigeo(i);
            xna[j] = fMesh->Reference()->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j);
            
            for (int k = DIM; k--; ){
                fIntegData.fAxes0(j,k) += xna[j] * dphigeo(k,i);
            };
        };
    };

    int ncols = DIM;
    int dim = DIM;

    switch (DIM) {
        case 2:
        {
            // //Computing the jacobian determinant and Inverse
            fIntegData.fA0(0,0) = fIntegData.fAxes0(0,0);
            fIntegData.fA0(0,1) = fIntegData.fAxes0(0,1);
            fIntegData.fA0(1,0) = fIntegData.fAxes0(1,0);
            fIntegData.fA0(1,1) = fIntegData.fAxes0(1,1);
            fIntegData.fJacA0 = fIntegData.fA0(0,0) * fIntegData.fA0(1,1) - fIntegData.fA0(0,1) * fIntegData.fA0(1,0);

            fIntegData.fA0Inv(0,0) = fIntegData.fA0(1,1) / fIntegData.fJacA0;
            fIntegData.fA0Inv(1,1) = fIntegData.fA0(0,0) / fIntegData.fJacA0;
            fIntegData.fA0Inv(0,1) = -fIntegData.fA0(0,1) / fIntegData.fJacA0;
            fIntegData.fA0Inv(1,0) = -fIntegData.fA0(1,0) / fIntegData.fJacA0;

            // std::cout << "fIntegData.fJacA0 = " << fIntegData.fJacA0 << std::endl;
            // std::cout << "fIntegData.fA0 = \n" << fIntegData.fA0 << std::endl;
            // std::cout << "fIntegData.fA0Inv = \n" << fIntegData.fA0Inv << std::endl;

            fIntegData.fJacA0 = fabs(fIntegData.fJacA0);
#ifdef DEBUG_BUILD
            if (fIntegData.fJacA0 < 1.e-6) PanicButton();
#endif
        }
            break;
        default:
            PanicButton();
            break;
    }

    return;
};

template<class geohape>
void GeoElementT<geoshape>::ComputeCurrentJacobian() {

    int DIM = geoshape::Dimension;
    fIntegData.fA1.resize(DIM,DIM);
    fIntegData.fA1.setZero();
    fIntegData.fX1.resize(3);
    fIntegData.fX1.setZero();
    fIntegData.fAxes1.resize(3,DIM);
    fIntegData.fAxes1.setZero();
    if (fIntegData.fAxes1Prev.size() != 0){
        fIntegData.fAxes1Prev.setZero();
        fIntegData.fA1Prev.setZero();
    }

    VecDouble yna(3);
    yna.setZero();
    
    fIntegData.fA1.setZero();
    for (int i = geoshape::NShape; i--; ){
        for (int j = fWeakForm->NState(); j--; ){
            // Approximate the integration space
            fIntegData.fX1[j] += (fMesh->Reference()->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) + fMesh->ConnectVec()[fGeoNodes[i]] -> GetSolution(j)) * fIntegData.fPhi(i);
            yna[j] = fMesh->Reference()->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) + fMesh->ConnectVec()[fGeoNodes[i]] -> GetSolution(j);
            double yprev = 0.;
            if (fIntegData.fAxes1Prev.size() != 0){
                yprev = fMesh->Reference()->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) + fMesh->ConnectVec()[fGeoNodes[i]] -> GetPreviousSolution(j);
            }


            for (int k = DIM; k--; ){
                fIntegData.fAxes1(j,k) += yna[j] * fIntegData.fDPhi(k,i);
                if (fIntegData.fAxes1Prev.size() != 0){
                    fIntegData.fAxes1Prev(j,k) += yprev * fIntegData.fDPhi(k,i);
                }
            };
        };
    };

    switch (compshape::Dimension)
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

        if (fIntegData.fA1Prev.size() != 0){
            fIntegData.fA1Prev(0,0) = fIntegData.fAxes1Prev(0,0);
            fIntegData.fA1Prev(0,1) = fIntegData.fAxes1Prev(0,1);
            fIntegData.fA1Prev(1,0) = fIntegData.fAxes1Prev(1,0);
            fIntegData.fA1Prev(1,1) = fIntegData.fAxes1Prev(1,1);
        }
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
template<class geohape>
void GeoElementT<geoshape>::ComputeSpatialDerivatives() {
    fIntegData.fPhi.setZero();
    fIntegData.fDPhi.setZero();

    VecInt orders(compshape::NSides);
    for (int i = compshape::NSides; i--; ) orders[i] = this->fMesh->ConnectVec()[fConnect[i]] -> GetOrder();

    //Shape functions
    compshape::Shape(fIntegData.fAdimCoord,fIntegData.fPhi,fIntegData.fDPhi,orders);
    // shapeQuad.ShapeHessian(xsi,ddphi);

    //Shape functions spatial first derivatives
    fIntegData.fDPhiX0 = fIntegData.fA0Inv.transpose() * fIntegData.fDPhi;

    return;
};

template<class geohape>
void GeoElementT<geoshape>::ComputeCurrentSpatialDerivatives() {
    VecInt orders(compshape::NSides);
    for (int i = compshape::NSides; i--; ) orders[i] = this->fMesh->ConnectVec()[fConnect[i]] -> GetOrder();

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
template<class geohape>
void GeoElementT<geoshape>::ComputeHighOrderSpatialDerivatives() {
    
    int DIM = compshape::Dimension;
        
    int nshape = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder());

    dDphi_dx.setZero();
    std::vector<MatrixDouble> ddphi(nshape,MatrixDouble(DIM,DIM));
    for (int i = 0; i < nshape; i++)
    {
        ddphi[i].setZero();
    }
    
    compshape::ShapeHessian(xsi,ddphi);
    
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

    for (int i = nshape; i--; ){
        xna_.setZero();
        for (int j = DIM; j--; ){
            // Approximate the integration space
            xna_[j] = fMesh->Reference()->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j);
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
        for (int i = nshape; i--; ){
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






#include "ShapePoint.h"
#include "ShapeOneDLin.h"
#include "ShapeOneDQua.h"
#include "ShapeOneDCub.h"
#include "ShapeQuadrilateralLin.h"
#include "ShapeQuadrilateralQua.h"
#include "ShapeTetrahedronLin.h"
#include "ShapeTetrahedronQua.h"
#include "ShapeTetrahedronCub.h"
#include "ShapeTriangleLin.h"
#include "ShapeTriangleQua.h"
#include "ShapeTriangleCub.h"
#include "ShapeHexahedron.h"

template class GeoElementT<ShapePoint>;
template class GeoElementT<ShapeOneDLin>;
template class GeoElementT<ShapeOneDQua>;
template class GeoElementT<ShapeOneDCub>;
template class GeoElementT<ShapeQuadrilateralLin>;
template class GeoElementT<ShapeQuadrilateralQua>;
template class GeoElementT<ShapeTetrahedronLin>;
template class GeoElementT<ShapeTetrahedronQua>;
template class GeoElementT<ShapeTetrahedronCub>;
template class GeoElementT<ShapeTriangleLin>;
template class GeoElementT<ShapeTriangleQua>;
template class GeoElementT<ShapeTriangleCub>;
template class GeoElementT<ShapeHexahedron>;