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
    
    MatrixDouble fTangentTensor(6,6);

    fTangentTensor.setZero();
    
    MatrixDouble NxN = fFlowVector.TensorProduct(fFlowVector);
    
    double Afactor = 2. * fShearModulus *(1.-3.*fShearModulus*data.fPlasticMultiplier/fVonMisesStress);
    double Bfactor = 6. * fShearModulus * fShearModulus * (data.fPlasticMultiplier/fVonMisesStress - 1./(3.*fShearModulus+fHardening));

    //Elastic operator 7.107
    if (fPlaneStress){
        PanicButton();
    } else {
        //Elastic Tensor
        // fTangentTensor = 2.*fShearModulus*fIdentity4Dev + fBulkModulus*fId2xId2;
        //Elastoplastic tensor
        fTangentTensor = Afactor*fIdentity4Dev + Bfactor*NxN + fBulkModulus*fId2xId2;  
    }
    
    
    // auto elast = fElasticModel->ConstitutiveMatrix();
    // std::cout << "Elastic - \n" << elast << std::endl;
    // std::cout << "Tangent - \n" << fTangentTensor << std::endl;

    if (fElasticModel->Dimension() == 2){
        VecInt order(3);
        order[0] = XX;
        order[1] = YY;
        order[2] = XY;
        MatrixDouble fTangent2D(3,3);
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                fTangent2D(i,j) = fTangentTensor(order[i],order[j]);
            }
        }
        fElasticModel->ConstitutiveMatrix() = fTangent2D;
    } else if (fElasticModel->Dimension() == 3) {
        fElasticModel->ConstitutiveMatrix() = fTangentTensor;
    } else {
        PanicButton();
    }
    
    fElasticModel->ComputeStiffness(index,data,Stiffness);
};
    
void VonMises::ComputeError(IntPointData &data, VecDouble &errors){
    fElasticModel->ComputeError(data,errors);
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
        fDeviatory = Stress.Deviatory();
        fVonMisesStress = sqrt(1.5*fDeviatory.DoubleContraction(fDeviatory));
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
    auto hydrostatic = Stress.Hydrostatic();

    fFlowVector = FlowVector(Stress);

    fDeviatory *= (1. - data.fPlasticMultiplier * 3. * fShearModulus / fVonMisesStress);
    Tensor epsilonUpdated(fDeviatory);
    epsilonUpdated *= 1./(2. * fShearModulus);
    Tensor Ident;
    double epslion_e_trial = data.fElasticStrain[index].Trace();
    Ident.Identity();
    Ident *= epslion_e_trial / 3.;
    epsilonUpdated += Ident;

    // Tensor epsilonUpdated2(data.fElasticStrain[index]);
    // double devnorm = fDeviatory.Norm();
    // epsilonUpdated2 += fDeviatory * (-data.fPlasticMultiplier * sqrt(1.5) / devnorm);
    data.fElasticStrain[index] = epsilonUpdated;

    
    Stress = fDeviatory + hydrostatic;

    data.fPlasticStrain[index] += data.fPlasticMultiplier;
}   