#include "LinearHardening.h"
#include "ElasticTruss.h"
#include "Elasticity2D.h"

LinearHardening::LinearHardening(WeakForm *elast, double hardModulus, double yield) : PlasticityModel(elast){
    fHardening = hardModulus;
    fYield = yield;

    ElasticTruss *truss = dynamic_cast<ElasticTruss* >(fElasticModel);
    Elasticity2D *mat2d = dynamic_cast<Elasticity2D* >(fElasticModel);
    if (truss) {
        fConstitutiveMatrix = truss->ConstitutiveMatrix();
        fYoungModulus = truss->YoungModulus();
    }
    if (mat2d) {
        fConstitutiveMatrix = mat2d->ConstitutiveMatrix();
        fYoungModulus = mat2d->YoungModulus();
    }
}

void LinearHardening::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    fElasticModel->ComputeStiffness(index,data,Stiffness);
    // Change stiffness;
    if (fPlasticStrain.rows()>0){
        Stiffness *= fHardening / (fYoungModulus+fHardening);
    }
};
    
void LinearHardening::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){
    fElasticModel->ComputeResidual(index,data,Rhs);
    if (fPlasticStrain.rows()>0){
        if (fRealDimension == 1){
            int nphi = data.fPhi.size();
            double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
            double elementLenght = 2.*data.fJacA0;
            double K = (fYoungModulus / elementLenght);

            MatrixDouble rotation(fDimension*nphi,fDimension*nphi);
            MatrixDouble matB(fDimension,fDimension*nphi);
            rotation.setZero();
            matB.setZero();
            double cosa = data.fAxes0(0,0) / data.fJacA0;
            double sina = data.fAxes0(1,0) / data.fJacA0;
            double check = sina*sina+cosa*cosa;
            for (int j = 0; j < nphi; j++){
                // for (int i = 0; i < fDimension; i++){
                    matB(0,fDimension*j) = data.fDPhiX0(j,0);
                // }
                rotation(2*j  ,2*j  ) = cosa;
                rotation(2*j+1,2*j  ) = sina;
                rotation(2*j  ,2*j+1) = -sina;
                rotation(2*j+1,2*j+1) = cosa;
            }

            VecDouble sol(2);// = data.fSol;
            sol[0] = -fPlasticStrain(0,0)*sina;
            sol[1] = +fPlasticStrain(0,0)*cosa;

            // std::cout << "rotation =\n"<< rotation << std::endl;
            Rhs += rotation * matB.transpose() * sol * WJ * elementLenght * K;
        } else if (fRealDimension == 2){
            int nphi = data.fPhi.size();

            double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
            MatrixDouble matB(3,2*nphi);
            matB.setZero();
            
            for (int j = 0; j < nphi; j++){
                matB(0,fDimension*j  ) = data.fDPhiX0(j,0);
                matB(1,fDimension*j+1) = data.fDPhiX0(j,1);
                matB(2,fDimension*j  ) = data.fDPhiX0(j,1);
                matB(2,fDimension*j+1) = data.fDPhiX0(j,0);
            }
            
            VecDouble auxplasticstrain;
            TensorToVoigt(fPlasticStrain,auxplasticstrain);
            
            Rhs -= matB.transpose() * fConstitutiveMatrix * auxplasticstrain * WJ;
        } else {
            PanicButton();
        }
    }
    
};
    
void LinearHardening::ComputeError(IntPointData &data, VecDouble &errors){
    fElasticModel->ComputeError(data,errors);

};

int LinearHardening::VariableIndex(const std::string &name) const{
    return fElasticModel->VariableIndex(name);

};

int LinearHardening::NSolutionVariables(int var) const {
    return fElasticModel->NSolutionVariables(var);
};

void LinearHardening::Solution(IntPointData &data, int var, VecDouble &Sol){
    fElasticModel->Solution(data,var,Sol);
};

void LinearHardening::ComputePlasticStrain(IntPointData &data, MatrixDouble &plasticstrain, MatrixDouble &totalstrain){
    //Check if the Integration Point is in the elastic region
    //1 - Compute integration point stress
    int var = fElasticModel->VariableIndex("Stress");
    int nsol = fElasticModel->NSolutionVariables(var);
    VecDouble Sol(nsol);
    MatrixDouble StressTensor(fRealDimension,fRealDimension);
    StressTensor.setZero();
    fElasticModel->Solution(data,var,Sol);
    //Compute the principal stress'
    VoigtToTensor(StressTensor,Sol);
    VecDouble PrincipalStress(fRealDimension);
    ComputePrincipalStress(StressTensor,PrincipalStress);


    if (fabs(PrincipalStress[0]) > fYield + plasticstrain.norm() * fHardening){
        VecDouble auxStress, auxTotalStrain;
        TensorToVoigt(StressTensor,auxStress);
        TensorToVoigt(totalstrain,auxTotalStrain);
        VecDouble DeltaStrain = fConstitutiveMatrix.inverse() * auxStress - auxTotalStrain;
        VecDouble DeltaPlasticStrain = (fYoungModulus/(fYoungModulus+fHardening)) * DeltaStrain;
        VecDouble DeltaSigma = (fYoungModulus*fHardening/(fYoungModulus+fHardening)) * DeltaStrain;
        MatrixDouble auxDStrain, auxDPlastStrain;
        VoigtToTensor(auxDStrain,DeltaStrain);
        VoigtToTensor(auxDPlastStrain,DeltaPlasticStrain);
        totalstrain += auxDStrain;
        plasticstrain += auxDPlastStrain;
        fTotalStrain = totalstrain;
        fPlasticStrain = plasticstrain;
    } else {
        VecDouble auxStress;
        TensorToVoigt(StressTensor,auxStress);
        VecDouble auxStrain = fConstitutiveMatrix.inverse() * auxStress;
        VoigtToTensor(totalstrain,auxStrain);
        
        if (fPlasticStrain.rows() > 0){
            fPlasticStrain.resize(0,0);
            fTotalStrain.resize(0,0);
        }
        
    }

}