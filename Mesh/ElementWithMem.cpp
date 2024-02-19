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
        this->fIntegData.fElasticStrain.resize(this->fIntRule.NPoints());
        for (int i = 0; i < this->fIntRule.NPoints(); i++){
            this->fIntegData.fElasticStrain[i].Zero();
            // this->fIntegData.fElasticStrain[i].setZero();
        }
        fElasticConstitutiveMatrix = fPlasticityModel->ElasticModel()->ConstitutiveMatrix();
        
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
    this->fIntegData.fDSolDxPrev.resize(this->fWeakForm->NState(), DIM);
    this->fIntegData.fNeedsSol = true;
    this->fIntegData.fSol.resize(this->fWeakForm->NState());
    this->fIntegData.fSolPrev.resize(this->fWeakForm->NState());

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
        
        //Assemble the trial stress tensor
        //the first step compute the initial elastic strain and store it
        Tensor ElasStress;
        VecDouble stressVoigt;
        int dim = fPlasticityModel->RealDimension();
        MatrixDouble fElasticStress(dim,dim);
        fElasticStress.setZero();
        int var = fPlasticityModel->ElasticModel()->VariableIndex("DeltaStrain");
        int nsol = fPlasticityModel->ElasticModel()->NSolutionVariables(var);
        VecDouble Sol(nsol);
        Sol.setZero();
        fPlasticityModel->ElasticModel()->Solution(this->fIntegData,var,Sol);
        this->fIntegData.fElasticStrain[index].fXX() += Sol[0];
        this->fIntegData.fElasticStrain[index].fYY() += Sol[1];
        this->fIntegData.fElasticStrain[index].fXY() += Sol[2];
        VecDouble elasticStrainTrial(3);
        elasticStrainTrial[0] = this->fIntegData.fElasticStrain[index].fXX();
        elasticStrainTrial[1] = this->fIntegData.fElasticStrain[index].fYY();
        elasticStrainTrial[2] = this->fIntegData.fElasticStrain[index].fXY();
        // if (this->fIntegData.fElasticStrain[index].norm() == 0){
        //     stressVoigt = fElasticConstitutiveMatrix * Sol;
        //     this->fIntegData.fElasticStrain[index] = Sol;
        // } else {
            stressVoigt = fElasticConstitutiveMatrix * elasticStrainTrial;
        // }
        // if (this->fIntegData.fElasticStrain[index].norm() == 0){
            
        // }
        fPlasticityModel->VoigtToTensor(fElasticStress,stressVoigt);
        ElasStress.SetData(fElasticStress);
        //Add sigma z component for plane strain problem
        if (dim == 2 && !fPlasticityModel->PlaneStress()){
            double E = fPlasticityModel->YoungModulus();
            double nu = fPlasticityModel->PoissonRatio();
            double k = E / ((1.+nu)*(1.-2.*nu));
            ElasStress.fZZ() = k * (nu * Sol[0] + nu * Sol[1]);
        }
        // std::cout << "Elastic tangent = \n" << fElasticConstitutiveMatrix << std::endl;

        //Check the Yield crieterion
        double YieldFunction = fPlasticityModel->YieldFunction(index,this->fIntegData,ElasStress);

        if (YieldFunction < 1.e-10){
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