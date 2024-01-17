#include "ElementWithMem.h"
#include "PlasticityModel.h"
#include "Elasticity2D.h"
#include "ElasticityPositional2D.h"
#include "ElasticTruss.h"
#include "PositionalTruss.h"

template<class tshape>
ElementWithMem<tshape>::ElementWithMem(int64_t index, VecInt &connect, CompMesh* mesh, WeakForm *wf) : ElementT<tshape>(index,connect,mesh,wf){
    fPlasticityModel = dynamic_cast<PlasticityModel *> (wf);
    if (fPlasticityModel){
        int nintpoints = this->fIntRule.NPoints();
        fPlasticStrain.resize(nintpoints);
        fTotalStrain.resize(nintpoints);
        int realdim = fPlasticityModel->RealDimension();
        for (int i = 0; i < nintpoints; i++){
            fPlasticStrain[i].resize(realdim,realdim);
            fPlasticStrain[i].setZero();
            fTotalStrain[i].resize(realdim,realdim);
            fTotalStrain[i].setZero();
        }
    } else {
        auto connect = this->getConnectivity();
        int nstate = wf->NState();
        for (int i = 0; i < connect.size(); i++){
            this->Mesh()->NodeVec()[connect[i]]->SetNStateVariables(nstate);
        }
    }

};

template<class tshape>
void ElementWithMem<tshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

    if (!this->fWeakForm) return;

    int DIM = tshape::Dimension;
    
    int index = 0;
    this->fIntegData.fAdimCoord.resize(DIM);
    this->fIntegData.fNeedsDSol = true;
    this->fIntegData.fDSolDx.resize(fPlasticityModel->NState(), DIM);
    this->fIntegData.fNeedsSol = true;
    this->fIntegData.fSol.resize(fPlasticityModel->NState());
    
    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (fPlasticityModel->ElasticModel());
    auto *truss = dynamic_cast<PositionalTruss *> (fPlasticityModel->ElasticModel());

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

        //Update the plastic strain
        fPlasticityModel->ComputePlasticStrain(this->fIntegData,fPlasticStrain[it],fTotalStrain[it]);

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
template<class tshape>
void ElementWithMem<tshape>::ComputeElContribution(std::vector<MatrixDouble> &jacobianNRMatrix, std::vector<VecDouble> &rhsVector){

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

        //Update the plastic strain
        
        //Computes the element diffusion/viscosity matrix
        this->fWeakForm->ComputeStiffness(index, this->fIntegData, jacobianNRMatrix);

        if (this->fIntegData.fNeedsSol) this->interpolateSolution();
        if (this->fIntegData.fNeedsDSol) this->interpolateSolDerivatives();

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

template class ElementWithMem<ShapePoint>;
template class ElementWithMem<ShapeOneDLin>;
template class ElementWithMem<ShapeOneDQua>;
template class ElementWithMem<ShapeOneDCub>;
template class ElementWithMem<ShapeTriangleLin>;
template class ElementWithMem<ShapeTriangleQua>;
template class ElementWithMem<ShapeTriangleCub>;
template class ElementWithMem<ShapeQuadrilateralLin>;
template class ElementWithMem<ShapeQuadrilateralQua>;
template class ElementWithMem<ShapeTetrahedronLin>;
template class ElementWithMem<ShapeTetrahedronQua>;
template class ElementWithMem<ShapeTetrahedronCub>;
template class ElementWithMem<ShapeHexahedron>;