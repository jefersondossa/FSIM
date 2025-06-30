#include "ElementTransient.h"
#include "Elasticity2D.h"
#include "ElasticityPositional2D.h"
#include "ElasticTruss.h"
#include "PositionalTruss.h"

template<class tshape>
ElementTransient<tshape>::ElementTransient(int64_t index, VecInt &connect, CompMesh* mesh, WeakForm *wf) : ElementT<tshape>(index,connect,mesh,wf){
    int DIM = tshape::Dimension;
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

template<class tshape>
void ElementTransient<tshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

    if (!this->fWeakForm) return;

    int DIM = tshape::Dimension;
    
    int index = 0;
    
    
    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (this->fWeakForm);
    auto *truss = dynamic_cast<PositionalTruss *> (this->fWeakForm);

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
        if (pos2d || truss){
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

    // std::cout << "Stiffness \n" << jacobianNRMatrix << '\n';
    // std::cout << "Rhs \n" << rhsVector << '\n';

    return;
};


//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<class tshape>
void ElementTransient<tshape>::ComputeElContribution(std::vector<MatrixDouble> &jacobianNRMatrix, std::vector<VecDouble> &rhsVector){
    PanicButton();
    if (!this->fWeakForm) return;

    int DIM = tshape::Dimension;
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

template class ElementTransient<ShapePoint>;
template class ElementTransient<ShapeOneDLin>;
template class ElementTransient<ShapeOneDQua>;
template class ElementTransient<ShapeOneDCub>;
template class ElementTransient<ShapeTriangleLin>;
template class ElementTransient<ShapeTriangleQua>;
template class ElementTransient<ShapeTriangleCub>;
template class ElementTransient<ShapeQuadrilateralLin>;
template class ElementTransient<ShapeQuadrilateralQua>;
template class ElementTransient<ShapeTetrahedronLin>;
template class ElementTransient<ShapeTetrahedronQua>;
template class ElementTransient<ShapeTetrahedronCub>;
template class ElementTransient<ShapeHexahedron>;