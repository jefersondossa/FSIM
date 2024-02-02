#include "Tresca.h"
#include "ElasticTruss.h"
#include "Elasticity2D.h"
#include "PositionalTruss.h"
#include "ElasticityPositional2D.h"

Tresca::Tresca(WeakForm *elast) : PlasticityModel(elast){

    // ElasticTruss *truss = dynamic_cast<ElasticTruss* >(fElasticModel);
    // Elasticity2D *mat2d = dynamic_cast<Elasticity2D* >(fElasticModel);
    // PositionalTruss *postruss = dynamic_cast<PositionalTruss* >(fElasticModel);
    // if (truss) {
    //     fConstitutiveMatrix = truss->ConstitutiveMatrix();
    //     fYoungModulus = truss->YoungModulus();
    //     fPoissonRatio = 0.;
    // }
    // if (mat2d) {
    //     fConstitutiveMatrix = mat2d->ConstitutiveMatrix();
    //     fYoungModulus = mat2d->YoungModulus();
    //     fPoissonRatio = mat2d->PoissonRatio();
    // }
    // if (postruss) {
    //     fConstitutiveMatrix = postruss->ConstitutiveMatrix();
    //     fYoungModulus = postruss->YoungModulus();
    //     fPoissonRatio = 0.;
    // }
}

void Tresca::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    fElasticModel->ComputeStiffness(index,data,Stiffness);
 
};
    
void Tresca::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){
    fElasticModel->ComputeResidual(index,data,Rhs);
    
};
    
void Tresca::ComputeError(IntPointData &data, VecDouble &errors){
    fElasticModel->ComputeError(data,errors);

};

int Tresca::VariableIndex(const std::string &name) const{
    return fElasticModel->VariableIndex(name);

};

int Tresca::NSolutionVariables(int var) const {
    return fElasticModel->NSolutionVariables(var);
};

void Tresca::Solution(IntPointData &data, int var, VecDouble &Sol){
    fElasticModel->Solution(data,var,Sol);
};

void Tresca::ComputePlasticStrain(IntPointData &data, MatrixDouble &plasticstrain, MatrixDouble &totalstrain){
    // //Check if the Integration Point is in the elastic region
    // //1 - Compute integration point stress
    // int var = fElasticModel->VariableIndex("Stress");
    // int nsol = fElasticModel->NSolutionVariables(var);
    // VecDouble Sol(nsol), auxPlasticStrain;
    // MatrixDouble StressTensor(fRealDimension,fRealDimension);
    // StressTensor.setZero();
    // fElasticModel->Solution(data,var,Sol);
    // //Compute the principal stress'
    
    // TensorToVoigt(plasticstrain,auxPlasticStrain);
    // Sol -= fConstitutiveMatrix * auxPlasticStrain;
    // VoigtToTensor(StressTensor,Sol);
    // VecDouble PrincipalStress(fRealDimension);
    // ComputePrincipalStress(StressTensor,PrincipalStress);

    // double C1 = pow(StressTensor(0,0) + StressTensor(1,1),2.);
    // double C2 = pow(StressTensor(0,0) - StressTensor(1,1),2.) + StressTensor(0,1)*StressTensor(0,1);
    // // double C3 = fYoungModulus / (3.*(1.-fElasticModel->Y))

    // double f = StressTensor.norm() - sqrt(2./3.) * fHardening;

    // double maxStress = PrincipalStress.maxCoeff();
    // double minStress = PrincipalStress.minCoeff();
    // double maxAbsStress = std::max(fabs(maxStress),fabs(minStress));

    // if (maxAbsStress > fabs(f)){
    //     VecDouble auxStress, auxTotalStrain;
    //     TensorToVoigt(StressTensor,auxStress);
    //     TensorToVoigt(totalstrain,auxTotalStrain);
    //     VecDouble DeltaStrain = fConstitutiveMatrix.inverse() * auxStress - auxTotalStrain;
    //     VecDouble DeltaPlasticStrain = (fYoungModulus/(fYoungModulus+fHardening)) * DeltaStrain;
    //     VecDouble DeltaSigma = (fYoungModulus*fHardening/(fYoungModulus+fHardening)) * DeltaStrain;
    //     MatrixDouble auxDStrain, auxDPlastStrain;
    //     VoigtToTensor(auxDStrain,DeltaStrain);
    //     VoigtToTensor(auxDPlastStrain,DeltaPlasticStrain);
    //     totalstrain += auxDStrain;
    //     plasticstrain += auxDPlastStrain;
    //     fTotalStrain = totalstrain;
    //     fPlasticStrain = plasticstrain;
    // } else {
    //     VecDouble auxStress;
    //     TensorToVoigt(StressTensor,auxStress);
    //     VecDouble auxStrain = fConstitutiveMatrix.inverse() * auxStress;
    //     VoigtToTensor(totalstrain,auxStrain);
        
    //     if (fPlasticStrain.rows() > 0){
    //         fPlasticStrain.resize(0,0);
    //         fTotalStrain.resize(0,0);
    //     }
        
    // }

}