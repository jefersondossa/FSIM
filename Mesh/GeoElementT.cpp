#include "GeoElementT.h"

template<class geoshape>
GeoElementT<geoshape>::GeoElementT() : GeoElement(){

};

template<class geoshape>
GeoElementT<geoshape>::GeoElementT(int64_t index, VecInt &geonodes, GeoMesh* mesh, int &material) : GeoElement(){
    fMesh = mesh;
    fGeoNodes.resize(geoshape::NShape);
    fIndex = index;
    for (int i = geoshape::NShape; i--; ) fGeoNodes[i] = geonodes[i];
    fMaterial = material;
    fNeighborElements.clear();
};



//------------------------------------------------------------------------------
//-------------------------SPATIAL TRANSFORM - JACOBIAN-------------------------
//------------------------------------------------------------------------------
template<class geoshape>
void GeoElementT<geoshape>::ComputeJacobian(IntPointData &data) {

    int DIM = geoshape::Dimension;
    data.fA0Inv.setZero();
    data.fAxes0.setZero();
    data.fA0.setZero();
    data.fX.setZero();
    VecDouble phigeo(geoshape::NShape);
    MatrixDouble dphigeo(geoshape::Dimension,geoshape::NShape);
    VecInt orders(geoshape::NSides);
    orders.fill(1);

    geoshape::Shape(data.fAdimCoord,phigeo,dphigeo,orders);

    data.fA0.setZero();
    VecDouble xna(3);
    for (int i = geoshape::NShape; i--; ){
        for (int j = 3; j--; ){
            // Approximate the integration space
            data.fX[j] += fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) * phigeo(i);
            // std::cout << "Coord " << i << " " << j << " = " << fMesh->NodeVec()[fConnect[i]] -> getCoordinateValue(j) << std::endl;
            xna[j] = fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j);
            
            for (int k = DIM; k--; ){
                data.fAxes0(j,k) += xna[j] * dphigeo(k,i);
            };
        };
    };

    int ncols = DIM;
    int dim = DIM;

    switch (DIM) {
        case 0:
            data.fJacA0 = 1.;
            break;
        case 1:
        {
            VecDouble v_1(3);
            v_1.setZero();

            for (int i = 0; i < 3; i++) {
                v_1[i] = data.fAxes0(i, 0);
            }

            double norm_v_1 = 0.;
            for (int i = 0; i < 3; i++) {
                norm_v_1 += v_1[i] * v_1[i];
            }

            norm_v_1 = sqrt(norm_v_1);
            data.fA0(0, 0) = norm_v_1;
            data.fJacA0 = norm_v_1;
            data.fA0Inv(0, 0) = 1.0 / data.fJacA0;

            data.fJacA0 = fabs(data.fJacA0);

        }
            break;
        case 2:
        {
            // auto *BC = dynamic_cast<L2Projection *> (fWeakForm);

            if (fMesh->Dimension() == 3) {
                VecDouble v_1(3), v_2(3);
                VecDouble v_1_til(3), v_2_til(3);
                v_1.setZero();
                v_2.setZero();
                v_1_til.setZero();
                v_2_til.setZero();

                for (int i = 0; i < 3; i++) {
                    v_1[i] = data.fAxes0(i, 0);
                    v_2[i] = data.fAxes0(i, 1);
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


                data.fA0(0, 0) = norm_v_1_til;
                data.fA0(0, 1) = v_1_dot_v_2 / norm_v_1_til;
                data.fA0(1, 1) = norm_v_2_til;

                data.fJacA0 = data.fA0(0, 0) * data.fA0(1, 1) - data.fA0(1, 0) * data.fA0(0, 1);

                data.fA0Inv(0, 0) = +data.fA0(1, 1) / data.fJacA0;
                data.fA0Inv(1, 1) = +data.fA0(0, 0) / data.fJacA0;
                data.fA0Inv(0, 1) = -data.fA0(0, 1) / data.fJacA0;
                data.fA0Inv(1, 0) = -data.fA0(1, 0) / data.fJacA0;

                data.fJacA0 = fabs(data.fJacA0);

                for (int i = 0; i < 3; i++) {
                    v_2_til[i] /= norm_v_2_til;
                    data.fAxes0(i,0) = v_1_til[i];
                    data.fAxes0(i,1) = v_2_til[i];
                }

                // std::cout << "data.fJacA0 = " << data.fJacA0 << std::endl;
                // std::cout << "data.fA0 = \n" << data.fA0 << std::endl;
                // std::cout << "data.fA0Inv = \n" << data.fA0Inv << std::endl;
            } else {
                //Computing the jacobian determinant and Inverse
                data.fA0(0,0) = data.fAxes0(0,0);
                data.fA0(0,1) = data.fAxes0(0,1);
                data.fA0(1,0) = data.fAxes0(1,0);
                data.fA0(1,1) = data.fAxes0(1,1);
                data.fJacA0 = data.fA0(0,0) * data.fA0(1,1) - data.fA0(0,1) * data.fA0(1,0);

                data.fA0Inv(0,0) = data.fA0(1,1) / data.fJacA0;
                data.fA0Inv(1,1) = data.fA0(0,0) / data.fJacA0;
                data.fA0Inv(0,1) = -data.fA0(0,1) / data.fJacA0;
                data.fA0Inv(1,0) = -data.fA0(1,0) / data.fJacA0;

                // std::cout << "data.fJacA0 = " << data.fJacA0 << std::endl;
                // std::cout << "data.fA0 = \n" << data.fA0 << std::endl;
                // std::cout << "data.fA0Inv = \n" << data.fA0Inv << std::endl;

                data.fJacA0 = fabs(data.fJacA0);
            }
            
#ifdef DEBUG_BUILD
            if (data.fJacA0 < 1.e-6) PanicButton();
#endif
        }
            break;
        case 3:
        {
            // axes.resize(dim, 3);

            for (int i = 0; i < 3; i++) {
                data.fA0(i, 0) = data.fAxes0(i, 0);
                data.fA0(i, 1) = data.fAxes0(i, 1);
                data.fA0(i, 2) = data.fAxes0(i, 2);
            }
            data.fJacA0 = 0.;
            data.fJacA0 -= data.fA0(0, 2) * data.fA0(1, 1) * data.fA0(2, 0); //- a02 a11 a20
            data.fJacA0 += data.fA0(0, 1) * data.fA0(1, 2) * data.fA0(2, 0); //+ a01 a12 a20
            data.fJacA0 += data.fA0(0, 2) * data.fA0(1, 0) * data.fA0(2, 1); //+ a02 a10 a21
            data.fJacA0 -= data.fA0(0, 0) * data.fA0(1, 2) * data.fA0(2, 1); //- a00 a12 a21
            data.fJacA0 -= data.fA0(0, 1) * data.fA0(1, 0) * data.fA0(2, 2); //- a01 a10 a22
            data.fJacA0 += data.fA0(0, 0) * data.fA0(1, 1) * data.fA0(2, 2); //+ a00 a11 a22

            data.fA0Inv(0, 0) = (-data.fA0(1, 2) * data.fA0(2, 1) + data.fA0(1, 1) * data.fA0(2, 2)) / data.fJacA0; //-a12 a21 + a11 a22
            data.fA0Inv(0, 1) = (data.fA0(0, 2) * data.fA0(2, 1) - data.fA0(0, 1) * data.fA0(2, 2)) / data.fJacA0; //a02 a21 - a01 a22
            data.fA0Inv(0, 2) = (-data.fA0(0, 2) * data.fA0(1, 1) + data.fA0(0, 1) * data.fA0(1, 2)) / data.fJacA0; //-a02 a11 + a01 a12
            data.fA0Inv(1, 0) = (data.fA0(1, 2) * data.fA0(2, 0) - data.fA0(1, 0) * data.fA0(2, 2)) / data.fJacA0; //a12 a20 - a10 a22
            data.fA0Inv(1, 1) = (-data.fA0(0, 2) * data.fA0(2, 0) + data.fA0(0, 0) * data.fA0(2, 2)) / data.fJacA0; //-a02 a20 + a00 a22
            data.fA0Inv(1, 2) = (data.fA0(0, 2) * data.fA0(1, 0) - data.fA0(0, 0) * data.fA0(1, 2)) / data.fJacA0; //a02 a10 - a00 a12
            data.fA0Inv(2, 0) = (-data.fA0(1, 1) * data.fA0(2, 0) + data.fA0(1, 0) * data.fA0(2, 1)) / data.fJacA0; //-a11 a20 + a10 a21
            data.fA0Inv(2, 1) = (data.fA0(0, 1) * data.fA0(2, 0) - data.fA0(0, 0) * data.fA0(2, 1)) / data.fJacA0; //a01 a20 - a00 a21
            data.fA0Inv(2, 2) = (-data.fA0(0, 1) * data.fA0(1, 0) + data.fA0(0, 0) * data.fA0(1, 1)) / data.fJacA0; //-a01 a10 + a00 a11

            data.fJacA0 = fabs(data.fJacA0);

            // axes.setZero();
            // axes(0, 0) = 1.0;
            // axes(1, 1) = 1.0;
            // axes(2, 2) = 1.0;
        }
            break;
    }

    return;
};


template<class geoshape>
void GeoElementT<geoshape>::ComputeJacobianSearch(IntPointData &data) {

    int DIM = geoshape::Dimension;
    data.fA0Inv.setZero();
    data.fAxes0.setZero();
    data.fA0.setZero();
    data.fX.setZero();
    data.fPhi.setZero();
    data.fDPhi.setZero();
    VecDouble phigeo(geoshape::NShape);
    MatrixDouble dphigeo(geoshape::Dimension,geoshape::NShape);
    VecInt orders(geoshape::NSides);
    orders.fill(1);

    geoshape::Shape(data.fAdimCoord,phigeo,dphigeo,orders);
    data.fPhi = phigeo;
    data.fDPhi = dphigeo;
    data.fA0.setZero();
    VecDouble xna(3);
    for (int i = geoshape::NShape; i--; ){
        for (int j = 3; j--; ){
            // Approximate the integration space
            data.fX[j] += fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) * phigeo(i);
            xna[j] = fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j);
            
            for (int k = DIM; k--; ){
                data.fAxes0(j,k) += xna[j] * dphigeo(k,i);
            };
        };
    };

    int ncols = DIM;
    int dim = DIM;

    switch (DIM) {
        case 1:
        {
            data.fA0.resize(2,2);
            data.fA0Inv.resize(2,2);
            VecDouble v_1(3),v_1Normal(3);
            v_1(0) = data.fAxes0(0, 0);
            v_1(1) = data.fAxes0(1, 0);
            v_1(2) = data.fAxes0(2, 0);
            double norm_v_1 = v_1.norm();
            v_1 /= norm_v_1;

            v_1Normal(0) = -v_1(1);
            v_1Normal(1) = v_1(0);
            v_1Normal(2) = 0.0;


            // //Computing the jacobian determinant and Inverse
            data.fA0(0,0) = data.fAxes0(0,0);
            data.fA0(1,0) = data.fAxes0(1,0);
            data.fA0(0,1) = v_1Normal(0);            
            data.fA0(1,1) = v_1Normal(1);
            data.fJacA0 = data.fA0(0,0) * data.fA0(1,1) - data.fA0(0,1) * data.fA0(1,0);

            data.fA0Inv(0,0) = data.fA0(1,1) / data.fJacA0;
            data.fA0Inv(1,1) = data.fA0(0,0) / data.fJacA0;
            data.fA0Inv(0,1) = -data.fA0(0,1) / data.fJacA0;
            data.fA0Inv(1,0) = -data.fA0(1,0) / data.fJacA0;

            // std::cout << "data.fJacA0 = " << data.fJacA0 << std::endl;
            // std::cout << "data.fA0 = \n" << data.fA0 << std::endl;
            // std::cout << "data.fA0Inv = \n" << data.fA0Inv << std::endl;

            data.fJacA0 = fabs(data.fJacA0);


            break;
        }
        case 2:
        {
            // //Computing the jacobian determinant and Inverse
            data.fA0(0,0) = data.fAxes0(0,0);
            data.fA0(0,1) = data.fAxes0(0,1);
            data.fA0(1,0) = data.fAxes0(1,0);
            data.fA0(1,1) = data.fAxes0(1,1);
            data.fJacA0 = data.fA0(0,0) * data.fA0(1,1) - data.fA0(0,1) * data.fA0(1,0);

            data.fA0Inv(0,0) = data.fA0(1,1) / data.fJacA0;
            data.fA0Inv(1,1) = data.fA0(0,0) / data.fJacA0;
            data.fA0Inv(0,1) = -data.fA0(0,1) / data.fJacA0;
            data.fA0Inv(1,0) = -data.fA0(1,0) / data.fJacA0;

            // std::cout << "data.fJacA0 = " << data.fJacA0 << std::endl;
            // std::cout << "data.fA0 = \n" << data.fA0 << std::endl;
            // std::cout << "data.fA0Inv = \n" << data.fA0Inv << std::endl;

            data.fJacA0 = fabs(data.fJacA0);
#ifdef DEBUG_BUILD
            if (data.fJacA0 < 1.e-6) PanicButton();
#endif
        }
            break;
        default:
            PanicButton();
            break;
    }

    return;
};

template<class geoshape>
void GeoElementT<geoshape>::ComputeCurrentJacobian(IntPointData &data, Element* compel) {

    int DIM = geoshape::Dimension;
    data.fA1.resize(DIM,DIM);
    data.fA1.setZero();
    data.fX1.resize(3);
    data.fX1.setZero();
    data.fAxes1.resize(3,DIM);
    data.fAxes1.setZero();
    if (data.fAxes1Prev.size() != 0){
        data.fAxes1Prev.setZero();
        data.fA1Prev.setZero();
    }

    auto *cmesh = compel->Mesh();
    auto *wf = compel->GetWeakForm();

    VecDouble yna(3);
    yna.setZero();
    
    data.fA1.setZero();
    for (int i = geoshape::NShape; i--; ){
        for (int j = wf->NState(); j--; ){
            // Approximate the integration space
            data.fX1[j] += (fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) + cmesh->ConnectVec()[fGeoNodes[i]] -> GetSolution(j)) * data.fPhi(i);
            yna[j] = fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) + cmesh->ConnectVec()[fGeoNodes[i]] -> GetSolution(j);
            double yprev = 0.;
            if (data.fAxes1Prev.size() != 0){
                yprev = fMesh->NodeVec()[fGeoNodes[i]] -> getCoordinateValue(j) + cmesh->ConnectVec()[fGeoNodes[i]] -> GetPreviousSolution(j);
            }


            for (int k = DIM; k--; ){
                data.fAxes1(j,k) += yna[j] * data.fDPhi(k,i);
                if (data.fAxes1Prev.size() != 0){
                    data.fAxes1Prev(j,k) += yprev * data.fDPhi(k,i);
                }
            };
        };
    };

    switch (geoshape::Dimension)
    {
    case 1:
    {
        VecDouble v_1(3);
        v_1.setZero();

        for (int i = 0; i < 3; i++) {
            v_1[i] = data.fAxes1(i, 0);
        }

        double norm_v_1 = 0.;
        for (int i = 0; i < 3; i++) {
            norm_v_1 += v_1[i] * v_1[i];
        }
        norm_v_1 = sqrt(norm_v_1);
        data.fA1(0, 0) = norm_v_1;
        data.fJacA1 = norm_v_1;
        data.fJacA1 = fabs(data.fJacA1);
    }
        break;
    case 2:
    {
        //Computing the jacobian determinant and Inverse
        data.fA1(0,0) = data.fAxes1(0,0);
        data.fA1(0,1) = data.fAxes1(0,1);
        data.fA1(1,0) = data.fAxes1(1,0);
        data.fA1(1,1) = data.fAxes1(1,1);
        data.fJacA1 = data.fA1(0,0) * data.fA1(1,1) - data.fA1(0,1) * data.fA1(1,0);
        data.fJacA1 = fabs(data.fJacA1);

        if (data.fA1Prev.size() != 0){
            data.fA1Prev(0,0) = data.fAxes1Prev(0,0);
            data.fA1Prev(0,1) = data.fAxes1Prev(0,1);
            data.fA1Prev(1,0) = data.fAxes1Prev(1,0);
            data.fA1Prev(1,1) = data.fAxes1Prev(1,1);
        }
    }
        break;
    case 3:
    {
        // axes.resize(dim, 3);

        for (int i = 0; i < 3; i++) {
            data.fA1(i, 0) = data.fAxes1(i, 0);
            data.fA1(i, 1) = data.fAxes1(i, 1);
            data.fA1(i, 2) = data.fAxes1(i, 2);
        }

        data.fJacA1 -= data.fA1(0, 2) * data.fA1(1, 1) * data.fA1(2, 0); //- a02 a11 a20
        data.fJacA1 += data.fA1(0, 1) * data.fA1(1, 2) * data.fA1(2, 0); //+ a01 a12 a20
        data.fJacA1 += data.fA1(0, 2) * data.fA1(1, 0) * data.fA1(2, 1); //+ a02 a10 a21
        data.fJacA1 -= data.fA1(0, 0) * data.fA1(1, 2) * data.fA1(2, 1); //- a00 a12 a21
        data.fJacA1 -= data.fA1(0, 1) * data.fA1(1, 0) * data.fA1(2, 2); //- a01 a10 a22
        data.fJacA1 += data.fA1(0, 0) * data.fA1(1, 1) * data.fA1(2, 2); //+ a00 a11 a22

        data.fJacA1 = fabs(data.fJacA1);

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
//----------------------SET ELEMENT INTERSECTION PARAMETERS---------------------
//------------------------------------------------------------------------------
template<class geoshape>
void GeoElementT<geoshape>::setIntersectionParameters(VecDouble &x, VecDouble &X) {
    xK.resize(2);
    XK.resize(2);
    xK[0] = x[0]; xK[1] = x[1]; 
    XK[0] = X[0]; XK[1] = X[1];
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