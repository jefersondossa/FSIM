#include "ElementTransient.h"
#include "Elasticity2D.h"
#include "ElasticityPositional2D.h"
#include "ElasticTruss.h"
#include "PositionalTruss.h"
#include "TransientPositionalFrame2D.h"

template<class geoshape, class compshape>
ElementTransient<geoshape,compshape>::ElementTransient(int64_t index, VecInt &connect, CompMesh* mesh, WeakForm *wf) : ElementT<geoshape,compshape>(index,connect,mesh,wf){
    int DIM = compshape::Dimension;
    this->fIntegData.fAdimCoord.resize(DIM);
    this->fIntegData.fAdimCoord.setZero();

    this->fIntegData.fNeedsDSol = true;
    this->fIntegData.fDSolDx.resize(this->fWeakForm->NState(), DIM);
    this->fIntegData.fDSolDx.setZero();

    this->fIntegData.fDSolDxPrev.resize(this->fWeakForm->NState(), DIM);
    this->fIntegData.fDSolDxPrev.setZero();

    this->fIntegData.fNeedsSol = true;
    this->fIntegData.fSol.resize(this->fWeakForm->NState());
    this->fIntegData.fSol.setZero();

    this->fIntegData.fSolPrev.resize(this->fWeakForm->NState());
    this->fIntegData.fSolPrev.setZero();
    
    this->fIntegData.fDSolDt.resize(this->fWeakForm->NState());
    this->fIntegData.fDSolDt.setZero();

    this->fIntegData.fDSolDDt.resize(this->fWeakForm->NState());
    this->fIntegData.fDSolDDt.setZero();
    this->fIntegData.fNeedsTimeDerivatives = true;
};

template<class geoshape, class compshape>
void ElementTransient<geoshape,compshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

    if (!this->fWeakForm) return;

    int DIM = compshape::Dimension;

    int index = 0;
    
    
    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (this->fWeakForm);
    auto *pos2dt = dynamic_cast<TransientPositionalFrame2D *> (this->fWeakForm);
    auto *truss = dynamic_cast<PositionalTruss *> (this->fWeakForm);

    if (pos2dt){
        this->fIntegData.fDSolDAdim.resize(this->fWeakForm->NState(), DIM);
        this->fIntegData.fNeedsDSolDAdim = true;
    }

    for(int it = 0; it < this->fIntRule.NPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) this->fIntegData.fAdimCoord[k] = this->fIntRule.PointList(index,k);

        //Returns the quadrature integration weight
        this->fIntegData.fWeight = this->fIntRule.WeightList(index);

        //Computes the jacobian matrix
        this->ComputeJacobian();

        //Computes spatial derivatives
        this->ComputeSpatialDerivatives();

        // Computes current spatial derivatives (only for position-based weak forms)
        if (pos2d || pos2dt || truss){
            this->ComputeCurrentJacobian();
            this->ComputeCurrentSpatialDerivatives();
        }

        if (this->fIntegData.fNeedsSol) this->interpolateSolution();
        if (this->fIntegData.fNeedsDSol) this->interpolateSolDerivatives();
        if (this->fIntegData.fNeedsTimeDerivatives) this->interpolateSolDTimeDerivatives();

        //Computes the element diffusion/viscosity matrix
        this->fWeakForm->ComputeStiffness(index, this->fIntegData, jacobianNRMatrix);

        //Computes the RHS vector
        this->fWeakForm->ComputeResidual(index, this->fIntegData, rhsVector); 

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
void ElementTransient<geoshape,compshape>::ComputeElContribution(std::vector<MatrixDouble> &jacobianNRMatrix, std::vector<VecDouble> &rhsVector){
    PanicButton();
    if (!this->fWeakForm) return;

    int DIM = compshape::Dimension;
    this->fIntegData.fA0Inv.resize(DIM,DIM);
    this->fIntegData.fAdimCoord.resize(DIM);

    int index = 0;

    for(int it = 0; it < this->fIntRule.NPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int k = 0; k < DIM; k++) this->fIntegData.fAdimCoord[k] = this->fIntRule.PointList(index,k);

        //Returns the quadrature integration weight
        this->fIntegData.fWeight = this->fIntRule.WeightList(index);

        //Computes the jacobian matrix
        this->ComputeJacobian();

        //Computes spatial derivatives
        this->ComputeSpatialDerivatives();
        
        //Computes the element diffusion/viscosity matrix
        this->fWeakForm->ComputeStiffness(index, this->fIntegData, jacobianNRMatrix);

        if (this->fIntegData.fNeedsSol) this->interpolateSolution();
        if (this->fIntegData.fNeedsDSol) this->interpolateSolDerivatives();
        if (this->fIntegData.fNeedsTimeDerivatives) this->interpolateSolDTimeDerivatives();
        
        //Computes the RHS vector
        this->fWeakForm->ComputeResidual(index, this->fIntegData, rhsVector); 

        index++;        
    };  
    // // std::cout << "\nStiffness Element " << this->Index() << "\n" << jacobianNRMatrix[1];
    // // std::cout << "\nrhsVector Element " << this->Index() << "\n" << rhsVector[1];
    // //Apply boundary conditions
    // ApplyBC(jacobianNRMatrix, rhsVector);

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
#include "HierarquicalOneD.h"

template class ElementTransient<ShapePoint,ShapePoint>;
template class ElementTransient<ShapeOneDLin,ShapeOneDLin>;
template class ElementTransient<ShapeOneDQua,ShapeOneDQua>;
template class ElementTransient<ShapeOneDCub,ShapeOneDCub>;
template class ElementTransient<ShapeOneDLin,HierarquicalOneD>;
template class ElementTransient<ShapeOneDQua,HierarquicalOneD>;
template class ElementTransient<ShapeOneDCub,HierarquicalOneD>;
template class ElementTransient<ShapeTriangleLin,ShapeTriangleLin>;
template class ElementTransient<ShapeTriangleQua,ShapeTriangleQua>;
template class ElementTransient<ShapeTriangleCub,ShapeTriangleCub>;
template class ElementTransient<ShapeQuadrilateralLin,ShapeQuadrilateralLin>;
template class ElementTransient<ShapeQuadrilateralQua,ShapeQuadrilateralQua>;
template class ElementTransient<ShapeTetrahedronLin,ShapeTetrahedronLin>;
template class ElementTransient<ShapeTetrahedronQua,ShapeTetrahedronQua>;
template class ElementTransient<ShapeTetrahedronCub,ShapeTetrahedronCub>;
template class ElementTransient<ShapeHexahedron,ShapeHexahedron>;