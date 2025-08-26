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
        this->fIntegData.fPlasticMultiplier.resize(this->fIntRule.NPoints());
        this->fIntegData.fPlasticMultiplier.setZero();
        this->fIntegData.fEffectivePlasticStrain.resize(this->fIntRule.NPoints());
        this->fIntegData.fEffectivePlasticStrain.setZero();
        this->fIntegData.fElasticStrain.resize(this->fIntRule.NPoints());
        this->fIntegData.fElasticStress.resize(this->fIntRule.NPoints());
        this->fIntegData.fElasticStrainIncrement.resize(this->fIntRule.NPoints());
        for (int i = 0; i < this->fIntRule.NPoints(); i++){
            this->fIntegData.fElasticStrain[i].Zero();
            this->fIntegData.fElasticStress[i].Zero();
            this->fIntegData.fElasticStrainIncrement[i].Zero();
            // this->fIntegData.fElasticStrain[i].setZero();
        }
        fElasticConstitutiveMatrix.resize(6,6);
        fElasticConstitutiveMatrix = 2.*fPlasticityModel->ShearModulus()*fPlasticityModel->fIdentity4Dev + fPlasticityModel->BulkModulus()*fPlasticityModel->fId2xId2;
        if (fPlasticityModel->PlaneStress()){
            fElasticConstitutiveMatrix.setZero();
            double alpha = (3.*fPlasticityModel->BulkModulus() - 2.*fPlasticityModel->ShearModulus()) / (3.*fPlasticityModel->BulkModulus() + 4.*fPlasticityModel->ShearModulus());
            fElasticConstitutiveMatrix(0,0) = fElasticConstitutiveMatrix(1,1) = 1. + alpha;
            fElasticConstitutiveMatrix(0,1) = fElasticConstitutiveMatrix(1,0) = alpha;
            fElasticConstitutiveMatrix(5,5) = 0.5;
            fElasticConstitutiveMatrix *= 2.*fPlasticityModel->ShearModulus();
        }

        auto felastPos = dynamic_cast<ElasticityPositional2D *> (fPlasticityModel->ElasticModel());
        if (felastPos){
            this->fIntegData.fAxes1Prev.resize(3,2);
            this->fIntegData.fAxes1Prev.setZero();
            this->fIntegData.fA1Prev.resize(2,2);
            this->fIntegData.fA1Prev.setZero();
        }
    } else {
        PanicButton();
    }
    int DIM = tshape::Dimension;
    this->fIntegData.fAdimCoord.resize(DIM);
    this->fIntegData.fNeedsDSol = true;
    this->fIntegData.fDSolDx.resize(this->fWeakForm->NState(), DIM);
    this->fIntegData.fDSolDxPrev.resize(this->fWeakForm->NState(), DIM);
    this->fIntegData.fNeedsSol = true;
    this->fIntegData.fSol.resize(this->fWeakForm->NState());
    this->fIntegData.fSolPrev.resize(this->fWeakForm->NState());
};

template<class tshape>
void ElementWithMem<tshape>::ComputeTrialStress(int &index,Tensor3D &ElasStress){
    
    auto fPlasticityModel = dynamic_cast<PlasticityModel *> (this->fWeakForm);
    
    //Assemble the trial stress tensor
    MatrixDouble fElasticStress(3,3);
    fElasticStress.setZero();
    int var = fPlasticityModel->ElasticModel()->VariableIndex("DeltaStrain");
    int nsol = fPlasticityModel->ElasticModel()->NSolutionVariables(var);
    VecDouble Sol(nsol);
    Sol.setZero();
    this->fIntegData.fIndex = index;
    fPlasticityModel->ElasticModel()->Solution(this->fIntegData,var,Sol);

    if (fPlasticityModel -> Dimension() == 2){
        this->fIntegData.fElasticStrain[index].fXX() += Sol[0];
        this->fIntegData.fElasticStrain[index].fYY() += Sol[1];
        this->fIntegData.fElasticStrain[index].fXY() += Sol[2];
        this->fIntegData.fElasticStrainIncrement[index].fXX() = Sol[0];
        this->fIntegData.fElasticStrainIncrement[index].fYY() = Sol[1];
        this->fIntegData.fElasticStrainIncrement[index].fXY() = Sol[2];
        if (fPlasticityModel->PlaneStress()){
            double poisson = fPlasticityModel->PoissonRatio();
            this->fIntegData.fElasticStrain[index].fZZ() -= poisson/(1.-poisson) * (Sol[0] + Sol[1]);
            this->fIntegData.fElasticStrainIncrement[index].fZZ() = -poisson/(1.-poisson) * (Sol[0] + Sol[1]);
        }   
    } else {
        this->fIntegData.fElasticStrain[index].fXX() += Sol[0];
        this->fIntegData.fElasticStrain[index].fYY() += Sol[1];
        this->fIntegData.fElasticStrain[index].fZZ() += Sol[2];
        this->fIntegData.fElasticStrain[index].fYZ() += Sol[3];
        this->fIntegData.fElasticStrain[index].fXZ() += Sol[4];
        this->fIntegData.fElasticStrain[index].fXY() += Sol[5];
        this->fIntegData.fElasticStrainIncrement[index].fXX() = Sol[0];
        this->fIntegData.fElasticStrainIncrement[index].fYY() = Sol[1];
        this->fIntegData.fElasticStrainIncrement[index].fZZ() = Sol[2];
        this->fIntegData.fElasticStrainIncrement[index].fYZ() = Sol[3];
        this->fIntegData.fElasticStrainIncrement[index].fXZ() = Sol[4];
        this->fIntegData.fElasticStrainIncrement[index].fXY() = Sol[5];
    }
    ElasStress = (this->fIntegData.fElasticStrain[index]).Multiply(fElasticConstitutiveMatrix);
    // ElasStress.fData = fElasticConstitutiveMatrix * this->fIntegData.fElasticStrain[index].fData;
    // int a = 0.;
    if (fPlasticityModel->PlaneStress()){
        ElasStress.fZZ() = 0.;
        ElasStress.fYZ() = 0.;
        ElasStress.fXZ() = 0.;
    }
}


template<class tshape>
void ElementWithMem<tshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix, VecDouble &rhsVector){

    if (!this->fWeakForm) return;
    auto fPlasticityModel = dynamic_cast<PlasticityModel *> (this->fWeakForm);
    
    int index = 0;
    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (fPlasticityModel->ElasticModel());
    auto *truss = dynamic_cast<PositionalTruss *> (fPlasticityModel->ElasticModel());
    int DIM = tshape::Dimension;
    
    for(int it = 0; it < this->fIntRule.NPoints(); it++){
        if (fPlasticityModel->Dimension() == 2){
            VecInt order(3);
            order[0] = XX;
            order[1] = YY;
            order[2] = XY;
            MatrixDouble fTangent2D(3,3);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    fTangent2D(i,j) = fElasticConstitutiveMatrix(order[i],order[j]);
            fPlasticityModel->ElasticModel()->ConstitutiveMatrix() = fTangent2D;
        } else if (fPlasticityModel->Dimension() == 3) {
            fPlasticityModel->ElasticModel()->ConstitutiveMatrix() = fElasticConstitutiveMatrix;
        } else {
            PanicButton();
        }
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
        
        // Tensor3D ElasStress;
        ComputeTrialStress(index,this->fIntegData.fElasticStress[index]);        

        //Check the Yield crieterion
        double YieldFunction = fPlasticityModel->YieldFunction(index,this->fIntegData,this->fIntegData.fElasticStress[index]);

        if (YieldFunction < 1.e-8){
            //Elastic step
            fPlasticityModel->ElasticModel()->ComputeStiffness(index, this->fIntegData, jacobianNRMatrix);
            fPlasticityModel->ComputeResidual(index, this->fIntegData, rhsVector, this->fIntegData.fElasticStress[index]); 
        } else {
            //Plastic step
            this->fIntegData.fYieldFunction[index] = YieldFunction;
            this->fIntegData.fPlasticMultiplier[index] = fPlasticityModel->PlasticMultiplier(index,this->fIntegData,this->fIntegData.fElasticStress[index]);
            fPlasticityModel->UpdateStateVariables(index,this->fIntegData,this->fIntegData.fElasticStress[index]);
            fPlasticityModel->ComputeTangentStiffness(index, this->fIntegData, jacobianNRMatrix,this->fIntegData.fElasticStress[index]);
            if (pos2d || truss){
                fPlasticityModel->ElasticModel()->ComputeResidual(index, this->fIntegData, rhsVector, this->fIntegData.fElasticStress[index]); 
            }else{
                fPlasticityModel->ComputeResidual(index, this->fIntegData, rhsVector, this->fIntegData.fElasticStress[index]); 
            }
        }
        index++;
    };  

    // std::cout << "Stiffness \n" << jacobianNRMatrix << std::endl;
    // std::cout << "Rhs \n" << rhsVector << std::endl;

    return;
};


template<class tshape>
void ElementWithMem<tshape>::ComputeElContribution(MatrixDouble &jacobianNRMatrix){

    if (!this->fWeakForm) return;
    auto fPlasticityModel = dynamic_cast<PlasticityModel *> (this->fWeakForm);
    
    int index = 0;
    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (fPlasticityModel->ElasticModel());
    auto *truss = dynamic_cast<PositionalTruss *> (fPlasticityModel->ElasticModel());
    int DIM = tshape::Dimension;
    
    for(int it = 0; it < this->fIntRule.NPoints(); it++){
        if (fPlasticityModel->Dimension() == 2){
            VecInt order(3);
            order[0] = XX;
            order[1] = YY;
            order[2] = XY;
            MatrixDouble fTangent2D(3,3);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    fTangent2D(i,j) = fElasticConstitutiveMatrix(order[i],order[j]);
            fPlasticityModel->ElasticModel()->ConstitutiveMatrix() = fTangent2D;
        } else if (fPlasticityModel->Dimension() == 3) {
            fPlasticityModel->ElasticModel()->ConstitutiveMatrix() = fElasticConstitutiveMatrix;
        } else {
            PanicButton();
        }
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

        // if (this->fIntegData.fNeedsSol) this->interpolateSolution();
        // if (this->fIntegData.fNeedsDSol) this->interpolateSolDerivatives();
        
        // Tensor3D ElasStress;
        // ComputeTrialStress(index,ElasStress);        

        // //Check the Yield crieterion
        // double YieldFunction = fPlasticityModel->YieldFunction(index,this->fIntegData,ElasStress);
    
        if (this->fIntegData.fYieldFunction[index] < 1.e-8){
            //Elastic step
            fPlasticityModel->ElasticModel()->ComputeStiffness(index, this->fIntegData, jacobianNRMatrix);
        } else {
            //Plastic step
            // this->fIntegData.fYieldFunction[index] = YieldFunction;
            // this->fIntegData.fPlasticMultiplier[index] = fPlasticityModel->PlasticMultiplier(index,this->fIntegData,ElasStress);
            // fPlasticityModel->UpdateStateVariables(index,this->fIntegData,ElasStress);
            fPlasticityModel->ComputeTangentStiffness(index, this->fIntegData, jacobianNRMatrix,this->fIntegData.fElasticStress[index]);
            // if (pos2d || truss){
            //     fPlasticityModel->ElasticModel()->ComputeResidual(index, this->fIntegData, rhsVector, ElasStress); 
            // }else{
            //     fPlasticityModel->ComputeResidual(index, this->fIntegData, rhsVector, ElasStress); 
            // }
        }
        index++;
    };  

    // std::cout << "Stiffness \n" << jacobianNRMatrix << std::endl;
    // std::cout << "Rhs \n" << rhsVector << std::endl;

    return;
};



template<class tshape>
void ElementWithMem<tshape>::ComputeElContribution(VecDouble &rhsVector){

    if (!this->fWeakForm) return;
    auto fPlasticityModel = dynamic_cast<PlasticityModel *> (this->fWeakForm);
    
    int index = 0;
    auto *pos2d = dynamic_cast<ElasticityPositional2D *> (fPlasticityModel->ElasticModel());
    auto *truss = dynamic_cast<PositionalTruss *> (fPlasticityModel->ElasticModel());
    int DIM = tshape::Dimension;
    
    for(int it = 0; it < this->fIntRule.NPoints(); it++){
        if (fPlasticityModel->Dimension() == 2){
            VecInt order(3);
            order[0] = XX;
            order[1] = YY;
            order[2] = XY;
            MatrixDouble fTangent2D(3,3);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    fTangent2D(i,j) = fElasticConstitutiveMatrix(order[i],order[j]);
            fPlasticityModel->ElasticModel()->ConstitutiveMatrix() = fTangent2D;
        } else if (fPlasticityModel->Dimension() == 3) {
            fPlasticityModel->ElasticModel()->ConstitutiveMatrix() = fElasticConstitutiveMatrix;
        } else {
            PanicButton();
        }
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
        
        // Tensor3D ElasStress;
        ComputeTrialStress(index,this->fIntegData.fElasticStress[index]);        

        //Check the Yield crieterion
        double YieldFunction = fPlasticityModel->YieldFunction(index,this->fIntegData,this->fIntegData.fElasticStress[index]);

        if (YieldFunction < 1.e-8){
            //Elastic step
            fPlasticityModel->ComputeResidual(index, this->fIntegData, rhsVector, this->fIntegData.fElasticStress[index]); 
        } else {
            //Plastic step
            this->fIntegData.fYieldFunction[index] = YieldFunction;
            this->fIntegData.fPlasticMultiplier[index] = fPlasticityModel->PlasticMultiplier(index,this->fIntegData,this->fIntegData.fElasticStress[index]);
            fPlasticityModel->UpdateStateVariables(index,this->fIntegData,this->fIntegData.fElasticStress[index]);
            if (pos2d || truss){
                fPlasticityModel->ElasticModel()->ComputeResidual(index, this->fIntegData, rhsVector, this->fIntegData.fElasticStress[index]); 
            }else{
                fPlasticityModel->ComputeResidual(index, this->fIntegData, rhsVector, this->fIntegData.fElasticStress[index]); 
            }
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