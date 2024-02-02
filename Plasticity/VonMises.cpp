#include "VonMises.h"
#include "ElasticTruss.h"
#include "Elasticity2D.h"
#include "PositionalTruss.h"
#include "ElasticityPositional2D.h"

VonMises::VonMises(WeakForm *elast) : PlasticityModel(elast){
    ElasticTruss *truss = dynamic_cast<ElasticTruss* >(fElasticModel);
    Elasticity2D *mat2d = dynamic_cast<Elasticity2D* >(fElasticModel);
    PositionalTruss *postruss = dynamic_cast<PositionalTruss* >(fElasticModel);
    if (truss) {
        // fConstitutiveMatrix = truss->ConstitutiveMatrix();
        fYoungModulus = truss->YoungModulus();
        fPoissonRatio = 0.;
    }
    if (mat2d) {
        // fConstitutiveMatrix = mat2d->ConstitutiveMatrix();
        fYoungModulus = mat2d->YoungModulus();
        fPoissonRatio = mat2d->PoissonRatio();
    }
    if (postruss) {
        // fConstitutiveMatrix = postruss->ConstitutiveMatrix();
        fYoungModulus = postruss->YoungModulus();
    }
}

void VonMises::ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor &Stress){
    double G, K, H;
    H = 0.;
    
    int nSComponents = fRealDimension == 2 ? 3 : 6;
    VecDouble identity2(nSComponents), deviatory(nSComponents);
    MatrixDouble identity4(nSComponents,nSComponents);
    MatrixDouble identity4d(nSComponents,nSComponents);
    MatrixDouble fTangentTensor(nSComponents,nSComponents);
    fTangentTensor.setZero();
    identity2.setZero();
    identity4.setIdentity();
    identity4d.setZero();
    deviatory.setZero();
    auto devAux = Stress.Deviatory();

    if (fRealDimension == 2){   
        K = fYoungModulus / (2. * (1.-fPoissonRatio));
        G = fYoungModulus / (2. * (1.+fPoissonRatio));
        identity2[0] = 1.;
        identity2[1] = 1.;
        identity4(2,2) = 0.5;
        
        deviatory[0] = devAux.fData[0]; 
        deviatory[1] = devAux.fData[1]; 
        deviatory[2] = devAux.fData[4]; 

    } else if (fRealDimension == 3){
        K = fYoungModulus / (3. * (1.-2.*fPoissonRatio));
        G = fYoungModulus / (2. * (1.+fPoissonRatio));
        PanicButton();
    } else {
        PanicButton();
    }
    
    for (int i = 0; i < nSComponents; i++){
        for (int j = 0; j < nSComponents; j++){
            identity4d(i,j) = identity4(i,j) - identity2[i] * identity2[j] / 3.; 
        }
    }

    data.fPlasticMultiplier[index] /= (3. * G);
    double gamma = data.fPlasticMultiplier[index];
    data.fPlasticStrain[index] += gamma;

    double vonMisesStress = Stress.DeviatoryNorm();
    double Q = sqrt(1.5) * vonMisesStress;
    double Qtrial = Q + 3.*G*gamma;

    double Afactor = 2. * G *(1.-3.*G*gamma/Qtrial);
    double Bfactor = 6. * G * G * (gamma/Qtrial - 1./(3.*G+H))/(vonMisesStress*vonMisesStress);

    for (int i = 0; i < nSComponents; i++){
        for (int j = 0; j < nSComponents; j++){
            fTangentTensor(i,j) = Afactor*identity4d(i,j) + Bfactor*deviatory[i]*deviatory[j] + K*identity2[i]*identity2[j];  
        }
    }
    fElasticModel->ConstitutiveMatrix() = fTangentTensor;
    fElasticModel->ComputeStiffness(index,data,Stiffness);

};
    
void VonMises::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){
    fElasticModel->ComputeResidual(index,data,Rhs);
    // if (fPlasticStrain.rows()>0){
    //     if (fRealDimension == 1){
    //         int nphi = data.fPhi.size();
    //         double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    //         double elementLenght = 2.*data.fJacA0;
    //         double K = (fYoungModulus / elementLenght);

    //         MatrixDouble rotation(fDimension*nphi,fDimension*nphi);
    //         MatrixDouble matB(fDimension,fDimension*nphi);
    //         rotation.setZero();
    //         matB.setZero();
    //         double cosa = data.fAxes0(0,0) / data.fJacA0;
    //         double sina = data.fAxes0(1,0) / data.fJacA0;
    //         double check = sina*sina+cosa*cosa;
    //         for (int j = 0; j < nphi; j++){
    //             // for (int i = 0; i < fDimension; i++){
    //                 matB(0,fDimension*j) = data.fDPhiX0(j,0);
    //             // }
    //             rotation(2*j  ,2*j  ) = cosa;
    //             rotation(2*j+1,2*j  ) = sina;
    //             rotation(2*j  ,2*j+1) = -sina;
    //             rotation(2*j+1,2*j+1) = cosa;
    //         }

    //         VecDouble sol(2);// = data.fSol;
    //         sol[0] = -fPlasticStrain(0,0)*sina;
    //         sol[1] = +fPlasticStrain(0,0)*cosa;

    //         // std::cout << "rotation =\n"<< rotation << std::endl;
    //         Rhs -= rotation * matB.transpose() * sol * WJ * elementLenght * K;
    //     } else if (fRealDimension == 2){
    //         int nphi = data.fPhi.size();

    //         double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    //         MatrixDouble matB(3,2*nphi);
    //         matB.setZero();
            
    //         for (int j = 0; j < nphi; j++){
    //             matB(0,fDimension*j  ) = data.fDPhiX0(j,0);
    //             matB(1,fDimension*j+1) = data.fDPhiX0(j,1);
    //             matB(2,fDimension*j  ) = data.fDPhiX0(j,1);
    //             matB(2,fDimension*j+1) = data.fDPhiX0(j,0);
    //         }
            
    //         VecDouble auxplasticstrain;
    //         TensorToVoigt(fPlasticStrain,auxplasticstrain);
            
    //         Rhs -= matB.transpose() * fConstitutiveMatrix * auxplasticstrain * WJ;
    //     } else {
    //         PanicButton();
    //     }
    // }
    
};
    
void VonMises::ComputeError(IntPointData &data, VecDouble &errors){
    fElasticModel->ComputeError(data,errors);

};

int VonMises::VariableIndex(const std::string &name) const{
    return fElasticModel->VariableIndex(name);

};

int VonMises::NSolutionVariables(int var) const {
    return fElasticModel->NSolutionVariables(var);
};

void VonMises::Solution(IntPointData &data, int var, VecDouble &Sol){
    fElasticModel->Solution(data,var,Sol);
};

void VonMises::ComputePlasticStrain(IntPointData &data, MatrixDouble &plasticstrain, MatrixDouble &totalstrain){
    //Check if the Integration Point is in the elastic region
    //1 - Compute integration point Elastic stress
    int var = fElasticModel->VariableIndex("Stress");
    int nsol = fElasticModel->NSolutionVariables(var);
    VecDouble Sol(nsol), auxPlasticStrain;
    MatrixDouble fElasticStress(fRealDimension,fRealDimension);
    fElasticStress.setZero();
    fElasticModel->Solution(data,var,Sol);
    VoigtToTensor(fElasticStress,Sol);

    Tensor elastTensor(fElasticStress);

    //Compute the principal stress'
    



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

double VonMises::YieldFunction(int &index, IntPointData &data, Tensor &Stress){
    //Eq (7.74)
    double plasticstrain = data.fPlasticStrain[index];
    double fYield;
    fUniaxialYield(plasticstrain,fYield);
    double YF = sqrt(3.*Stress.J2()) - fYield;
    return YF;
}

Tensor VonMises::FlowVector(Tensor &Stress){
    //Eq (7.77)
    double devnorm = Stress.DeviatoryNorm();
    auto temp = Stress.Deviatory();
    temp *= sqrt(1.5) / devnorm;
    return temp;
}