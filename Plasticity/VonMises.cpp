#include "VonMises.h"

VonMises::VonMises(WeakForm *elast) : PlasticityModel(elast){

    if (fPlaneStress) {
        fMatP.resize(3,3);//Box 9.3
        fMatP.setZero();
        fMatP(0,0) = fMatP(1,1) = 2./3.;
        fMatP(0,1) = fMatP(1,0) =-1./3.;
        fMatP(2,2) = 6./3.;
        std::cout << "Please implement this option" << std::endl;
        PanicButton();
    }
}

void VonMises::ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor &Stress){
    
    VecDouble flowcurr(fNStressComponents);
    VecDouble flowPrev(fNStressComponents);
    MatrixDouble fTangentTensor(fNStressComponents,fNStressComponents);
    MatrixDouble fTangentTensor2(fNStressComponents,fNStressComponents);

    fTangentTensor.setZero();
    fTangentTensor2.setZero();
    auto devAux = Stress.Deviatory();
    auto flow = FlowVector(Stress);
    if (fRealDimension == 2){
        flowPrev[0] = fFlowVector.fXX(); 
        flowPrev[1] = fFlowVector.fYY(); 
        flowPrev[2] = fFlowVector.fXY();
        flowcurr[0] = flow.fXX(); 
        flowcurr[1] = flow.fYY(); 
        flowcurr[2] = flow.fXY();
    } else if (fRealDimension == 3){
        PanicButton();
    } else {
        PanicButton();
    }
    // 
    double q = Stress.DeviatoryNorm() ;//* sqrt(1.5);
    double qtrial = Stress.DeviatoryNorm() * sqrt(1.5) + 3.*fShearModulus*data.fPlasticMultiplier;

    double Afactor = 2. * fShearModulus *(1.-3.*fShearModulus*data.fPlasticMultiplier/fVonMisesStress);
    double Bfactor = 6. * fShearModulus * fShearModulus * (data.fPlasticMultiplier/fVonMisesStress - 1./(3.*fShearModulus+fHardening));

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            if (fPlaneStress){

            } else {
                fTangentTensor(i,j) = Afactor*fIdentity4Dev(i,j) + Bfactor*flowPrev[i]*flowPrev[j] + fBulkModulus*fIdentity2[i]*fIdentity2[j];  
                fTangentTensor2(i,j) = 2.*fShearModulus*fIdentity4Dev(i,j) + fBulkModulus*fIdentity2[i]*fIdentity2[j];// - 6. * fShearModulus * fShearModulus/ (3.*fShearModulus+ fHardening) * flowcurr[i]*flowcurr[j];
                // fTangentTensor2(i,j) = fConstitutiveMatrix(i,j) - 6. * fShearModulus * fShearModulus/ (3.*fShearModulus+ fHardening) * flow.fData[i]*flow.fData[j];
                // fTangentTensor(i,j) = 2.*fShearModulus*fIdentity4Dev(i,j) + fBulkModulus*fIdentity2[i]*fIdentity2[j] - 6.*fShearModulus*fShearModulus*flowcurr[i]*flowcurr[j]/(3.*fShearModulus+fHardening)
                //                     - data.fPlasticMultiplier * 6. * fShearModulus * fShearModulus/fVonMisesStress * (fIdentity4Dev(i,j) - flowPrev[i]*flowPrev[j]);  
            }
        }
    }

    //Elastic operator 7.107
    // auto elast = fElasticModel->ConstitutiveMatrix();
    // std::cout << "Elastic - \n" << elast << std::endl;
    // std::cout << "Tangent - \n" << fTangentTensor << std::endl;
    // std::cout << "Tangent2 - \n" << fTangentTensor2 << std::endl;

    // std::cout << "Differenca = \n" << elast - fTangentTensor2 << std::endl; 

    fElasticModel->ConstitutiveMatrix() = fTangentTensor;
    fElasticModel->ComputeStiffness(index,data,Stiffness);

};
    
void VonMises::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs, Tensor &Stress){

    // fElasticModel->ComputeResidual(index,data,Rhs);
    // return;
    int nphi = data.fPhi.size();

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    MatrixDouble matB(3,2*nphi);
    matB.setZero();
    for (int j = 0; j < nphi; j++){
        matB(0,2*j  ) = data.fDPhiX0(j,0);
        matB(1,2*j+1) = data.fDPhiX0(j,1);
        matB(2,2*j  ) = data.fDPhiX0(j,1);
        matB(2,2*j+1) = data.fDPhiX0(j,0);
    }

    auto force = fForceFunction;
    VecDouble forcingF(fDimension);
    forcingF.setZero();
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);
    
    VecDouble stress(3);

    // TensorToVoigt(Stress,stress);
    stress[0] = Stress.fXX();
    stress[1] = Stress.fYY();
    stress[2] = Stress.fXY();

    // Rhs -= matB.transpose() * fElasticModel->ConstitutiveMatrix() * strain * WJ;
    Rhs -= matB.transpose() * stress * WJ;

    for (int i = nphi; i--; ){
        double shapeFi = data.fPhi[i];
        //External force
        double Fx = forcingF[0] * shapeFi;
        double Fy = forcingF[1] * shapeFi;
        Rhs[2*i  ] += Fx * WJ;
        Rhs[2*i+1] += Fy * WJ;
    };
    
};
    
void VonMises::ComputeError(IntPointData &data, VecDouble &errors){
    fElasticModel->ComputeError(data,errors);

};

int VonMises::VariableIndex(const std::string &name) const{
    auto n = fElasticModel->VariableIndex(name);
    if (n != -1){
        return n;
    } else {
        if(!strcmp("PlasticStrain",name.c_str()))    return 101;
    }
    return -1;
    // return fElasticModel->VariableIndex(name);
};

int VonMises::NSolutionVariables(int var) const {
    int val = fElasticModel->NSolutionVariables(var);
    if (val != -1){
        return val;
    } else {
        switch (var)
        {
        case 101:
            return 1;
        
        default:
            return -1;
        }
    }
    return -1;
};

void VonMises::Solution(IntPointData &data, int var, VecDouble &Sol){
    fElasticModel->Solution(data,var,Sol);
    if (var == 101){
        if (data.fPlasticStrain.norm() > 0)
        Sol[0] = data.fPlasticStrain.mean();
        return;
    };
};
double VonMises::YieldFunction(int &index, IntPointData &data, Tensor &Stress){
    double YF = 0.;
    if (fPlaneStress){
        //Box 9.3
        VecDouble Sigma(3);
        Sigma[0] = Stress.fData[0];
        Sigma[1] = Stress.fData[1];
        Sigma[2] = Stress.fData[3];
        double sigmay;
        fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
        YF = 0.5 * Sigma.transpose()*fMatP*Sigma-sigmay*sigmay/3.; 
    } else {
        //Eq (7.74)
        double sigmay;
        fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
        auto deviatory = Stress.Deviatory();
        YF = sqrt(3.*Stress.J2()) - sigmay;
    }
    return YF;
}

Tensor VonMises::FlowVector(Tensor &Stress){
    //Eq (7.77)
    double devnorm = Stress.DeviatoryNorm();
    auto temp = Stress.Deviatory();
    temp *= 1. / devnorm;
    auto norm = temp.Norm();
    return temp;
}

double VonMises::PlasticMultiplier(int &index, IntPointData &data, Tensor &Stress){
    //Newton-Raphson to find plastic multiplier
    double dGamma = 0.;
    if (fPlaneStress){//Box 9.5
        VecDouble principalStress = Stress.Eigenvalues();
        double xi = (Stress.fXX()+Stress.fYY())*(Stress.fXX()+Stress.fYY()) / 6.
                  + (Stress.fYY()-Stress.fXX())*(Stress.fYY()-Stress.fXX()) / 2.
                  + 2. * Stress.fXY() * Stress.fXY();
        double sigmay;
        fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
        double PhiTil = 0.5 * xi - sigmay*sigmay / 3.;

    } else { //Box 7.4 Souza Neto
        auto deviatory = Stress.Deviatory();
        fVonMisesStress = sqrt(1.5*deviatory.DoubleContraction(deviatory));
        double sigmay = 0.;
        fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
        double deltaGamma = 0.;
        double PhiTil = fVonMisesStress-sigmay;
        
        while (fabs(PhiTil/sigmay) > 1.e-5){
            double d = -3. * fShearModulus - fHardening;
            dGamma -= PhiTil/d;
            fUniaxialYield(data.fPlasticStrain[index]+dGamma,sigmay,fHardening);
            PhiTil = fVonMisesStress - 3.*fShearModulus*dGamma-sigmay;
        }
    }
    double dgamma2 = data.fYieldFunction[index]/(3.*fShearModulus + fHardening);
    return dGamma;
}

void VonMises::UpdateStateVariables(int &index, IntPointData &data, Tensor &Stress){
    //Box 7.3
    auto dev = Stress.Deviatory();
    auto hydrostatic = Stress.Hydrostatic();
    fFlowVector = FlowVector(Stress);
    // implementar atualização da deformação elasticaç
    // Tensor devUpdated(dev);
    dev *= (1. - data.fPlasticMultiplier * 3. * fShearModulus / fVonMisesStress);
    Tensor epsilonUpdated(dev);
    epsilonUpdated *= 1./(2. * fShearModulus);
    Tensor Ident;
    double epslion_e_trial = data.fElasticStrain[index].Trace();
    Ident.Identity();
    Ident *= epslion_e_trial / 3.;
    epsilonUpdated += Ident;

    if (fNStressComponents == 2){
        data.fElasticStrain[index] = epsilonUpdated;
        // data.fElasticStrain[index][0] = epsilonUpdated.fXX();
        // data.fElasticStrain[index][1] = epsilonUpdated.fYY();
        // data.fElasticStrain[index][2] = epsilonUpdated.fXY();
    }
    
    // dev *= (1. - data.fPlasticMultiplier * 3 * fShearModulus / fVonMisesStress);
    Stress = dev + hydrostatic;
    data.fPlasticStrain[index] += data.fPlasticMultiplier;
}   