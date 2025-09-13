#include "ElementT.h"
#include "PositionalTruss.h"
#include "ElasticityPositional2D.h"
#include "L2Projection.h"
#include "PositionalFrame2D.h"

template<class geoshape, class compshape>
ElementT<geoshape,compshape>::ElementT() : Element(){
};

template<class geoshape, class compshape>
ElementT<geoshape,compshape>::ElementT(int64_t index, VecInt &geonodes, CompMesh* mesh, WeakForm *wf) : Element(){
    int DIM = compshape::Dimension;
    fMesh = mesh;
    fGeoNodes.resize(geoshape::NSides);
    fConnect.resize(compshape::NSides);
    fIndex = index;
    for (int i = geoshape::NSides; i--; ) fGeoNodes[i] = geonodes[i];
    DEG = fMesh->GetDefaultOrder();
    fWeakForm = wf;
    if (fWeakForm) nLocDOF = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()) * fWeakForm->NState();

    fNeighborElements.clear();
    int increase = 0;
    if(wf->GetExactSolution()) increase = 2;

    fIntRule.SetOrder(2*this->fMesh->GetDefaultOrder()+increase);
    // fIntRule.SetOrder(2);

    fIntegData.fWeightFunction.resize(fIntRule.NPoints());
    fIntegData.fDistFunction.resize(fIntRule.NPoints());
    fIntegData.fPrevWeightFunction.resize(fIntRule.NPoints());

    fIntegData.fWeightFunction.fill(1.);
    fIntegData.fPrevWeightFunction.fill(1.);

    getIntegPointCoordinates();
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

    fIntegData.fDPhi.resize(geoshape::Dimension,compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()));
    fIntegData.fDPhi.setZero();

    //Checks if it is a BC element and set the nodes as BC
    L2Projection *bc = dynamic_cast<L2Projection*>(fWeakForm);
    if (bc){
        if (bc->GetBCType() == BoundaryConditionType::kDirectionalHomogeneousDirichlet || 
            bc->GetBCType() == BoundaryConditionType::kDirichlet) {
        for (int i = 0; i < compshape::NSides; i++){
            fMesh->NodeVec()[fGeoNodes[i]]->SetHasBC();
        }
            }
    }
};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::setIntegPointWeightFunction() {
    //It produces wrong results for constant weight functions
    VecDouble xsi(geoshape::Dimension);    

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
           fIntegData.fWeightFunction[index] += fIntegData.fPhi[j] * fMesh->NodeVec()[fGeoNodes[j]] -> getWeightFunction();
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


template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeIntPointDistFunction(VecDouble &nodalval) {
    
    int DIM = geoshape::Dimension;
    fIntegData.fAdimCoord.resize(DIM);
    fIntegData.fPhi.resize(geoshape::NShape);
    
    // for(int i = 0; i < nQuad.getNumberOfIntegrationPoints(); i++) {
    //     intPointWeightFunctionPrev[i] = intPointWeightFunction[i];
    //     intPointWeightFunction[i] = 0.;
    // }

    int index=0;
    fIntegData.fDistFunction.setZero();
    for(int it = 0; it < fIntRule.NPoints(); it++){
        
        for (int i=0; i<DIM; i++) fIntegData.fAdimCoord[i] = fIntRule.PointList(index,i);
           
        compshape::Shape(fIntegData.fAdimCoord,fIntegData.fPhi);

        for (int j=0; j<geoshape::NShape; j++){
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
template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::getIntegPointCoordinates(){
    int DIM = compshape::Dimension;
    DEG = fMesh->GetDefaultOrder();

    fIntPointCoordinates.resize(fIntRule.NPoints(),2);
    fIntPointCoordinates.setZero();

    int nshape = geoshape::NShape;
    VecDouble xsi(DIM);
    VecDouble phi_(nshape);
    fIntPointCoordinates.resize(fIntRule.NPoints(),DIM);

    for (int i = 0; i < fIntRule.NPoints(); i++){
        double x[DIM] = {};

        for (int k = DIM; k--; ) xsi[k] = fIntRule.PointList(i,k);

        compshape::Shape(xsi,phi_);

        for (int k = DIM; k--; ) fIntPointCoordinates(i,k) = 0.;

        for (int j = 0; j < nshape; j++)
            for (int k = DIM; k--; )
                fIntPointCoordinates(i,k) += fMesh->NodeVec()[fGeoNodes[j]] -> getCoordinateValue(k) * phi_[j];
        
    };

    return;
};


//------------------------------------------------------------------------------
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------
template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::setIntersectionParameters(VecDouble &x, VecDouble &X) {
    xK.resize(2);
    XK.resize(2);
    xK[0] = x[0]; xK[1] = x[1]; 
    XK[0] = X[0]; XK[1] = X[1];
    return;
};


//------------------------------------------------------------------------------
//-------------------------SPATIAL TRANSFORM - JACOBIAN-------------------------
//------------------------------------------------------------------------------
template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeJacobian() {

    int DIM = geoshape::Dimension;
    fIntegData.fA0Inv.setZero();
    fIntegData.fAxes0.setZero();
    fIntegData.fA0.setZero();
    fIntegData.fX.setZero();
    fIntegData.fPhi.setZero();
    fIntegData.fDPhi.setZero();

    geoshape::Shape(fIntegData.fAdimCoord,fIntegData.fPhi);
    geoshape::ShapeGradient(fIntegData.fAdimCoord,fIntegData.fDPhi);   
   
    fIntegData.fA0.setZero();
    VecDouble xna(3);
    for (int i = geoshape::NShape; i--; ){
        for (int j = 3; j--; ){
            // Approximate the integration space
            fIntegData.fX[j] += fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) * fIntegData.fPhi(i);
            // std::cout << "Coord " << i << " " << j << " = " << fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j) << std::endl;
            xna[j] = fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j);
            
            for (int k = DIM; k--; ){
                fIntegData.fAxes0(j,k) += xna[j] * fIntegData.fDPhi(k,i);
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


template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeJacobianSearch() {

    int DIM = geoshape::Dimension;
    fIntegData.fA0Inv.setZero();
    fIntegData.fAxes0.setZero();
    fIntegData.fA0.setZero();
    fIntegData.fX.setZero();
    fIntegData.fPhi.setZero();
    fIntegData.fDPhi.setZero();

    geoshape::Shape(fIntegData.fAdimCoord,fIntegData.fPhi);
    geoshape::ShapeGradient(fIntegData.fAdimCoord,fIntegData.fDPhi);   

    fIntegData.fA0.setZero();
    VecDouble xna(3);
    for (int i = geoshape::NShape; i--; ){
        for (int j = 3; j--; ){
            // Approximate the integration space
            fIntegData.fX[j] += fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) * fIntegData.fPhi(i);
            xna[j] = fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j);
            
            for (int k = DIM; k--; ){
                fIntegData.fAxes0(j,k) += xna[j] * fIntegData.fDPhi(k,i);
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

template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeCurrentJacobian() {

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
            fIntegData.fX1[j] += (fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) + fMesh->ConnectVec()[fGeoNodes[i]] -> GetSolution(j)) * fIntegData.fPhi(i);
            yna[j] = fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) + fMesh->ConnectVec()[fGeoNodes[i]] -> GetSolution(j);
            double yprev = 0.;
            if (fIntegData.fAxes1Prev.size() != 0){
                yprev = fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) + fMesh->ConnectVec()[fGeoNodes[i]] -> GetPreviousSolution(j);
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
template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeSpatialDerivatives() {
    
    MatrixDouble DphiComp(compshape::Dimension,compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()));

    compshape::ShapeGradient(fIntegData.fAdimCoord,DphiComp);
    // shapeQuad.ShapeHessian(xsi,ddphi);

    //Shape functions spatial first derivatives
    fIntegData.fDPhiX0 = fIntegData.fA0Inv.transpose() * DphiComp;

    return;
};

template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeCurrentSpatialDerivatives() {
    
    
    MatrixDouble DphiComp(compshape::Dimension,compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()));

    compshape::ShapeGradient(fIntegData.fAdimCoord,DphiComp);
    // shapeQuad.ShapeHessian(xsi,ddphi);
    
    fIntegData.fDPhiX1.setZero();

    //Shape functions spatial first derivatives
    fIntegData.fDPhiX1 = fIntegData.fA1.inverse().transpose() * DphiComp;

    return;
};
//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::getHighOrderSpatialDerivatives(VecDouble &xsi, MatrixDouble &ainv_, MatrixDouble &dphi_dx, MatrixDouble &dDphi_dx) {
    
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
            xna_[j] = fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j);
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


//------------------------------------------------------------------------------
//---------------INTERPOLATES MESH VELOCITY AND ITS DERIVATIVES-----------------
//------------------------------------------------------------------------------
template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::interpolateMeshVelocity(int &index, VecDouble &umesh_, VecDouble &umeshPrev_) {

    int DIM = compshape::Dimension;
    umesh_.setZero();
    umeshPrev_.setZero();
    
    for (int i = geoshape::NShape; i--; ){
        double shapeFi = fIntegData.fPhi[i];
        for (int j = DIM; j--; ){
            umesh_[j] += fMesh->NodeVec()[fGeoNodes[i]] -> getMeshVelocity(j) * shapeFi;
        }
    }

    return;
}

template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::interpolateSolution(int &index, VecDouble &u_) {
    u_.setZero();

    for (int i = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()); i--; ){
        double shapeFi = fIntegData.fPhi[i];
        int nstate = fMesh->ConnectVec()[fConnect[i]]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            u_[j] += fMesh->ConnectVec()[fConnect[i]] -> GetSolution(j) * shapeFi;
        }
    }
}
template<class geoshape, class compshape>
double ElementT<geoshape,compshape>::InterpolateVariable(VecDouble &nValues, int point) {
    double val = 0.;
    int nshape = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder());
    fIntegData.fPhi.resize(nshape);
    fIntegData.fPhi.setZero();
    fIntegData.fAdimCoord.resize(compshape::Dimension);
    for (int i=0; i<compshape::Dimension; i++) fIntegData.fAdimCoord[i] = fIntRule.PointList(point,i);
    compshape::Shape(fIntegData.fAdimCoord,fIntegData.fPhi);
    for (int i = nshape; i--; ){
        double shapeFi = fIntegData.fPhi[i];
        val += nValues[i] * shapeFi;
    }
    return val;
}

template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::interpolateSolDTimeDerivatives() {
    fIntegData.fDSolDt.setZero();
    fIntegData.fDSolDDt.setZero();
    for (int i = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()); i--; ){
        double shapeFi = fIntegData.fPhi[i];
        int nstate = fMesh->ConnectVec()[fConnect[i]]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            fIntegData.fDSolDt[j] += fMesh->ConnectVec()[fConnect[i]] -> GetDSolutionDTime(j) * shapeFi;
            fIntegData.fDSolDDt[j] += fMesh->ConnectVec()[fConnect[i]] -> GetDSolutionDDTime(j) * shapeFi;
        }
    }
}

template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::interpolateSolDTimeDerivatives(VecDouble &du_dt, VecDouble &du_ddt) {
    fIntegData.fDSolDt.setZero();
    fIntegData.fDSolDDt.setZero();
    for (int i = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()); i--; ){
        double shapeFi = fIntegData.fPhi[i];
        int nstate = fMesh->ConnectVec()[fConnect[i]]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            fIntegData.fDSolDt[j] += fMesh->ConnectVec()[fConnect[i]] -> GetDSolutionDTime(j) * shapeFi;
            fIntegData.fDSolDDt[j] += fMesh->ConnectVec()[fConnect[i]] -> GetDSolutionDDTime(j) * shapeFi;
        }
    }
    du_dt = fIntegData.fDSolDt;
    du_ddt = fIntegData.fDSolDDt;
}


template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::interpolateSolution() {
    fIntegData.fSol.setZero();
    if (fIntegData.fSolPrev.size() != 0) fIntegData.fSolPrev.setZero();
    for (int i = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()); i--; ){
        double shapeFi = fIntegData.fPhi[i];
        int nstate = fMesh->ConnectVec()[fConnect[i]]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            fIntegData.fSol[j] += fMesh->ConnectVec()[fConnect[i]] -> GetSolution(j) * shapeFi;
            if (fIntegData.fSolPrev.size() != 0) {
                fIntegData.fSolPrev[j] += fMesh->ConnectVec()[fConnect[i]] -> GetPreviousSolution(j) * shapeFi;
            }
        }
    }

    // Store solution for nodes only in case the weakForm has no memory (used by topology optimization only) 
    if (fWeakForm && !fWeakForm->GetHasMemory()) {
        const size_t n_state = fMesh->ConnectVec()[fConnect[0]]->GetNStateVariables();

        fIntegData.fSolNodes.resize(geoshape::NShape * n_state);
        fIntegData.fSolNodes.setZero();
        for (int i = 0; i < geoshape::NShape; i++){
            for (int j = 0; j < n_state; j++ ){
                fIntegData.fSolNodes[(n_state*i) + j] = fMesh->ConnectVec()[fConnect[i]]->GetSolution(j);
            }
        }
    }
}

template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::interpolateSolution(VecDouble &phi, VecDouble &u_) {
    u_.setZero();
    for (int i = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()); i--; ){
        double shapeFi = phi(i);
        int nstate = fMesh->ConnectVec()[fConnect[i]]->GetNStateVariables();
        for (int j = 0; j < nstate; j++ ){
            u_[j] += fMesh->ConnectVec()[fConnect[i]] -> GetSolution(j) * shapeFi;
        }
    }
}

//------------------------------------------------------------------------------
//----------------------------INTERPOLATES VELOCITY-----------------------------
//------------------------------------------------------------------------------
template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::interpolateSolDerivatives(MatrixDouble &du_dx) {
    du_dx.setZero();    
    int DIM = compshape::Dimension;
    for (int i = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()); i--; ){
        int nstate = fMesh->ConnectVec()[fConnect[i]]->GetNStateVariables();
        for (int j = DIM; j--; ){
            for (int k = nstate; k--; ){
                du_dx(k,j) += fMesh->ConnectVec()[fConnect[i]] -> GetSolution(k) * fIntegData.fDPhiX0(j,i);
            }
        }
    }
}

template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::interpolateSolDerivatives(MatrixDouble &dphidx, MatrixDouble &du_dx) {
    du_dx.setZero();    
    int DIM = compshape::Dimension;
    for (int i = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()); i--; ){
        int nstate = fMesh->ConnectVec()[fConnect[i]]->GetNStateVariables();
        for (int j = DIM; j--; ){
            for (int k = nstate; k--; ){
                du_dx(k,j) += fMesh->ConnectVec()[fConnect[i]] -> GetSolution(k) * dphidx(j,i);
            }
        }
    }
}

template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::interpolateSolDerivatives() {
    fIntegData.fDSolDx.setZero();
    fIntegData.fDSolDAdim.setZero();
    bool flag = false;
    if (fIntegData.fDSolDxPrev.size() != 0){
        fIntegData.fDSolDxPrev.setZero(); 
        flag = true;
    }     
    for (int i = compshape::NShapeFunctions(this->fMesh->GetDefaultOrder()); i--; ){
        int nstate = fMesh->ConnectVec()[fConnect[i]]->GetNStateVariables();
        for (int j = compshape::Dimension; j--; ){
            for (int k = nstate; k--; ){
                fIntegData.fDSolDx(k,j) += fMesh->ConnectVec()[fConnect[i]] -> GetSolution(k) * fIntegData.fDPhiX0(j,i);
                if (flag) fIntegData.fDSolDxPrev(k,j) += fMesh->ConnectVec()[fConnect[i]] -> GetPreviousSolution(k) * fIntegData.fDPhiX0(j,i);

                if (fIntegData.fNeedsDSolDAdim) fIntegData.fDSolDAdim(k,j) += fMesh->ConnectVec()[fConnect[i]] -> GetSolution(k) * fIntegData.fDPhi(j,i);
#ifdef DEBUG_BUILD
                if (std::isnan(fIntegData.fDSolDx(k,j))){
                    PanicButton();
                }
#endif
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
template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

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
        ComputeJacobian();

        //Computes spatial derivatives
        ComputeSpatialDerivatives();
        
        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || truss || frame2d){
            ComputeCurrentJacobian();
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

    // std::cout << std::scientific << "Stiffness \n" << jacobianNRMatrix << std::endl;
    // std::cout << std::scientific << "Rhs \n" << rhsVector << std::endl;
    // Set stiffness matrix to cache.
    fIntegData.fStiffnessMatrix = jacobianNRMatrix;
    fIntegData.fRHS = rhsVector;

    return;
};

template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix){

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
        
        index++;        
    };  

    // std::cout << "Stiffness \n" << jacobianNRMatrix << std::endl;
    // std::cout << "Rhs \n" << rhsVector << std::endl;

    return;
};

template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeElContribution(VecDouble &rhsVector){

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
template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeElContribution(std::vector<MatrixDouble> &jacobianNRMatrix, std::vector<VecDouble> &rhsVector){

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
template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeElContribution(std::vector<MatrixDouble> &jacobianNRMatrix){

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
template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeElContribution(std::vector<VecDouble> &rhsVector){

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
        ComputeJacobian();

        //Computes spatial derivatives
        ComputeSpatialDerivatives();
        
        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || truss){
            ComputeCurrentJacobian();
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
template<class geoshape, class compshape>
void ElementT<geoshape,compshape>::ComputeError(VecDouble &errors){

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

template<class geoshape, class compshape >
Element * ElementT<geoshape,compshape>::Clone() const {
    return new ElementT(*this);
}

template<class geoshape, class compshape >
void ElementT<geoshape,compshape>::Integrate(std::vector<std::string> &varNames, std::map<std::string,VecDouble> &result){
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
            ComputeJacobian();

            //Computes spatial derivatives
            ComputeSpatialDerivatives();
            
            // Computes current spatial derivatives (only for position-based weak forms)
            // if (pos2d || truss){
            //     ComputeCurrentJacobian();
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

template class ElementT<ShapePoint,ShapePoint>;
template class ElementT<ShapeOneDLin,ShapeOneDLin>;
template class ElementT<ShapeOneDQua,ShapeOneDQua>;
template class ElementT<ShapeOneDCub,ShapeOneDCub>;
template class ElementT<ShapeOneDLin,HierarquicalOneD>;
template class ElementT<ShapeOneDQua,HierarquicalOneD>;
template class ElementT<ShapeOneDCub,HierarquicalOneD>;
template class ElementT<ShapeTriangleLin,ShapeTriangleLin>;
template class ElementT<ShapeTriangleQua,ShapeTriangleQua>;
template class ElementT<ShapeTriangleCub,ShapeTriangleCub>;
template class ElementT<ShapeQuadrilateralLin,ShapeQuadrilateralLin>;
template class ElementT<ShapeQuadrilateralQua,ShapeQuadrilateralQua>;
template class ElementT<ShapeTetrahedronLin,ShapeTetrahedronLin>;
template class ElementT<ShapeTetrahedronQua,ShapeTetrahedronQua>;
template class ElementT<ShapeTetrahedronCub,ShapeTetrahedronCub>;
template class ElementT<ShapeHexahedron,ShapeHexahedron>;