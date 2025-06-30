#include "ElementLagrangeMultiplier.h"
#include "Elasticity2D.h"
#include "ElasticityPositional2D.h"
#include "ElasticTruss.h"
#include "PositionalTruss.h"
#include "LagrangeMultiplier.h"

template<class tshape>
ElementLagrangeMultiplier<tshape>::ElementLagrangeMultiplier(int64_t index, Element* leftEl, Element* rightEl, CompMesh* cmesh, WeakForm *wf) : ElementT<tshape>(){
    this->fIndex = index;
    fLeftElement = leftEl;
    fRightElement = rightEl;
    this->fWeakForm = wf;
};

template<class tshape>
void ElementLagrangeMultiplier<tshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

    if (!this->fWeakForm) return;

    LagrangeMultiplier *lagrange = dynamic_cast<LagrangeMultiplier*>(this->fWeakForm);
    if (!lagrange) {
        PanicButton();
        return;
    }

    int DIM = tshape::Dimension;
    int index = 0;
    // auto intrule = fLeftElement->GetIntRule();
    
    for(int it = 0; it < fLeftElement->getNumberOfIntegrationPoints(); it++){

        // //Defines the integration points adimentional coordinates
        VecDouble Coord = fLeftElement->GetIntPointCoord(index);
        for (int k = 0; k < DIM; k++) fLeftElement->IntegrationData().fAdimCoord[k] = Coord[k];
        if (fRightElement->Dimension() != fLeftElement->Dimension()){
            fRightElement->IntegrationData().fAdimCoord.setZero();
        } else {
            for (int k = 0; k < DIM; k++) fRightElement->IntegrationData().fAdimCoord[k] = Coord[k];
        }
        

        // //Returns the quadrature integration weight
        fLeftElement->IntegrationData().fWeight = fLeftElement->GetIntPointWeight(index);
        fRightElement->IntegrationData().fWeight = fLeftElement->GetIntPointWeight(index);

        //Computes the jacobian matrix
        fLeftElement->ComputeJacobian();
        fRightElement->ComputeJacobian();

        //Computes spatial derivatives
        fLeftElement->ComputeSpatialDerivatives();
        fRightElement->ComputeSpatialDerivatives();

        if (fLeftElement->IntegrationData().fNeedsSol) fLeftElement->interpolateSolution();
        if (fLeftElement->IntegrationData().fNeedsDSol) fLeftElement->interpolateSolDerivatives();
        if (fRightElement->IntegrationData().fNeedsSol) fRightElement->interpolateSolution();
        if (fRightElement->IntegrationData().fNeedsDSol) fRightElement->interpolateSolDerivatives();

        //Computes the element diffusion/viscosity matrix
        lagrange->ComputeStiffness(index, fLeftElement->IntegrationData(), fRightElement->IntegrationData(), jacobianNRMatrix);

        //Computes the RHS vector
        lagrange->ComputeResidual(index, fLeftElement->IntegrationData(), fRightElement->IntegrationData(), rhsVector); 

        index++;        
    };  

    // std::cout << "Stiffness \n" << jacobianNRMatrix << '\n';
    // std::cout << "Rhs \n" << rhsVector << '\n';

    return;
};



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

template class ElementLagrangeMultiplier<ShapePoint>;
template class ElementLagrangeMultiplier<ShapeOneDLin>;
template class ElementLagrangeMultiplier<ShapeOneDQua>;
template class ElementLagrangeMultiplier<ShapeOneDCub>;
template class ElementLagrangeMultiplier<ShapeTriangleLin>;
template class ElementLagrangeMultiplier<ShapeTriangleQua>;
template class ElementLagrangeMultiplier<ShapeTriangleCub>;
template class ElementLagrangeMultiplier<ShapeQuadrilateralLin>;
template class ElementLagrangeMultiplier<ShapeQuadrilateralQua>;
template class ElementLagrangeMultiplier<ShapeTetrahedronLin>;
template class ElementLagrangeMultiplier<ShapeTetrahedronQua>;
template class ElementLagrangeMultiplier<ShapeTetrahedronCub>;
template class ElementLagrangeMultiplier<ShapeHexahedron>;