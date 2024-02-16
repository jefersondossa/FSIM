#include "ElementWithMem.h"
#include "PlasticityModel.h"
#include "Elasticity2D.h"
#include "ElasticityPositional2D.h"
#include "ElasticTruss.h"
#include "PositionalTruss.h"

template<class tshape>
ElementWithMem<tshape>::ElementWithMem(int64_t index, VecInt &connect, CompMesh* mesh, WeakForm *wf) : ElementT<tshape>(index,connect,mesh,wf){
    auto fPlasticityModel = dynamic_cast<PlasticityModel *> (wf);
    if (fPlasticityModel){
        this->fIntegData.fYieldFunction.resize(this->fIntRule.NPoints());
        this->fIntegData.fYieldFunction.setZero();
        this->fIntegData.fPlasticStrain.resize(this->fIntRule.NPoints());
        this->fIntegData.fPlasticStrain.setZero();
        this->fIntegData.fElasticStrain.resize(fPlasticityModel->NStressComponents());
        this->fIntegData.fElasticStrain.setZero();
    } else {

    }

};

template<class tshape>
void ElementWithMem<tshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

    if (!this->fWeakForm) return;
    auto fPlasticityModel = dynamic_cast<PlasticityModel *> (this->fWeakForm);

    int DIM = tshape::Dimension;
    
    int index = 0;
    this->fIntegData.fAdimCoord.resize(DIM);
    this->fIntegData.fNeedsDSol = true;
    this->fIntegData.fDSolDx.resize(this->fWeakForm->NState(), DIM);
    this->fIntegData.fNeedsSol = true;
    this->fIntegData.fSol.resize(this->fWeakForm->NState());
    
    
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

        //Assemble the stress tensor
        int var = fPlasticityModel->ElasticModel()->VariableIndex("Stress");
        int nsol = fPlasticityModel->ElasticModel()->NSolutionVariables(var);
        VecDouble Sol(nsol);
        Sol.setZero();
        int dim = fPlasticityModel->RealDimension();
        MatrixDouble fElasticStress(dim,dim);
        fElasticStress.setZero();
        fPlasticityModel->ElasticModel()->Solution(this->fIntegData,var,Sol);
        fPlasticityModel->VoigtToTensor(fElasticStress,Sol);

        
        Tensor ElasStress(fElasticStress);
        //Add sigma z component for plane strain problem
        if (dim == 2 && !fPlasticityModel->PlaneStress()){
            int varz = fPlasticityModel->ElasticModel()->VariableIndex("SigmaZ");
            int nsolz = fPlasticityModel->ElasticModel()->NSolutionVariables(varz);
            VecDouble Solz(nsolz);
            Solz.setZero();
            fPlasticityModel->ElasticModel()->Solution(this->fIntegData,varz,Solz);
            ElasStress.fZZ() = Solz[0];
        }

        //Check the Yield crieterion
        double YieldFunction = fPlasticityModel->YieldFunction(index,this->fIntegData,ElasStress);

        if (YieldFunction < 0){
            //Elastic step
            fPlasticityModel->ElasticModel()->ComputeStiffness(index, this->fIntegData, jacobianNRMatrix);
            fPlasticityModel->ElasticModel()->ComputeResidual(index, this->fIntegData, rhsVector); 
        } else {
            //Plastic step
            this->fIntegData.fYieldFunction[index] = YieldFunction;
            this->fIntegData.fPlasticMultiplier = fPlasticityModel->PlasticMultiplier(index,this->fIntegData,ElasStress);
            fPlasticityModel->UpdateStateVariables(index,this->fIntegData,ElasStress);
            fPlasticityModel->ComputeTangentStiffness(index, this->fIntegData, jacobianNRMatrix,ElasStress);
            fPlasticityModel->ComputeResidual(index, this->fIntegData, rhsVector, ElasStress); 
        }
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