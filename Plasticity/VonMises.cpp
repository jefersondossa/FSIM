#include "VonMises.h"

VonMises::VonMises(WeakForm *elast) : PlasticityModel(elast){
    fFlowVector.Zero();
    fDeviatory.Zero();

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
    
    double Afactor = 2. * fShearModulus *(1.-3.*fShearModulus*data.fPlasticMultiplier[index]/fVonMisesStress);
    double Bfactor = 6. * fShearModulus * fShearModulus * (data.fPlasticMultiplier[index]/fVonMisesStress - 1./(3.*fShearModulus+fHardening));

    //Elastic operator 7.107
    if (fPlaneStress){
        PanicButton();
    } else {
        //Elastic Tensor
        fTangentTensor = 2.*fShearModulus*fIdentity4Dev + fBulkModulus*fId2xId2;
        //Elastoplastic tensor
        //Continuum operator
        // fTangentTensor = 2.*fShearModulus*fIdentity4Dev + fBulkModulus*fId2xId2 - 6.*fShearModulus*fShearModulus/(3*fShearModulus+fHardening)*NxN;
        //Consistent operator
        // fTangentTensor = Afactor*fIdentity4Dev + Bfactor*NxN + fBulkModulus*fId2xId2;  
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
        if (fTangent2D.determinant() < 1.e-10){
            std::cout << "Elastoplastic operator with problem " << std::endl;
            PanicButton();
        }
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
        YF = sqrt(3.*Stress.J2()) - sigmay;
    }
    return YF;
}

Tensor VonMises::FlowVector(Tensor &Stress){
    //Eq (7.117)
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
        fVonMisesStress = sqrt(1.5*fDeviatory.DoubleContraction(fDeviatory))+3.*data.fPlasticMultiplier[index];
        double sigmay = 0.;
        fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
        double PhiTil = fVonMisesStress-sigmay;
        
        while (fabs(PhiTil/sigmay) > 1.e-5){
            double d = -3. * fShearModulus - fHardening;
            dGamma -= PhiTil/d;
            fUniaxialYield(data.fPlasticStrain[index]+dGamma,sigmay,fHardening);
            PhiTil = fVonMisesStress - 3.*fShearModulus*dGamma-sigmay;
        }
    }
    // dGamma = data.fYieldFunction[index]/(3.*fShearModulus + fHardening);
    //Update plastic strain
    data.fPlasticStrain[index] += dGamma;
    return dGamma;
}

void VonMises::UpdateStateVariables(int &index, IntPointData &data, Tensor &Stress){
    
    // auto elasticTrial = data.fElasticStrain[index];//+data.fElasticStrainIncrement[index];

    // double elasticStrainV = elasticTrial.Trace();
    // double pressure = fBulkModulus * elasticStrainV;

    // auto elasticTrialDev = elasticTrial.Deviatory();
    // elasticTrial.fXY() /= 2.;
    // elasticTrial.fXZ() /= 2.;
    // elasticTrial.fYZ() /= 2.;

    // double varJ2t=4.*fShearModulus*fShearModulus*(elasticTrialDev.fXY()*elasticTrialDev.fXY()+
    // 0.5*(elasticTrialDev.fXX()*elasticTrialDev.fXX()+elasticTrialDev.fYY()*elasticTrialDev.fYY()+elasticTrialDev.fZZ()*elasticTrialDev.fZZ()));

    // fFlowVector = FlowVector(elasticTrialDev);

    // double qtrial = sqrt(3.*varJ2t);
    // double sigmay;
    // fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);

    // //Newton-Raphson to find plastic multiplier
    // double dGamma = 0.;
    // if (fPlaneStress){//Box 9.5
    //     VecDouble principalStress = Stress.Eigenvalues();
    //     double xi = (Stress.fXX()+Stress.fYY())*(Stress.fXX()+Stress.fYY()) / 6.
    //               + (Stress.fYY()-Stress.fXX())*(Stress.fYY()-Stress.fXX()) / 2.
    //               + 2. * Stress.fXY() * Stress.fXY();
    //     double sigmay;
    //     fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
    //     double PhiTil = 0.5 * xi - sigmay*sigmay / 3.;

    // } else { //Box 7.4 Souza Neto
    //     fDeviatory = Stress.Deviatory();
    //     fVonMisesStress = sqrt(1.5*fDeviatory.DoubleContraction(fDeviatory))+3.*data.fPlasticMultiplier[index];
    //     double sigmay = 0.;
    //     fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
    //     double PhiTil = fVonMisesStress-sigmay;
    //     double ddGamma = 0.;
    //     while (fabs(PhiTil/sigmay) > 1.e-5){
    //         double d = -3. * fShearModulus - fHardening;
    //         ddGamma = - PhiTil/d;
    //         dGamma += ddGamma;
    //         data.fPlasticStrain[index] += ddGamma;
    //         fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
    //         PhiTil = fVonMisesStress - 3.*fShearModulus*dGamma-sigmay;
    //     }
    // }
    // data.fPlasticMultiplier[index] = dGamma;

    // double factor = 2.*fShearModulus*(1.-3.*fShearModulus*dGamma/qtrial);
    // Stress.fXX() = factor * elasticTrialDev.fXX() + pressure;
    // Stress.fYY() = factor * elasticTrialDev.fYY() + pressure;
    // Stress.fZZ() = factor * elasticTrialDev.fZZ() + pressure;
    // Stress.fXY() = factor * elasticTrialDev.fXY();

    // factor = factor/(2.*fShearModulus);

    // data.fElasticStrain[index].fXX() = factor*elasticTrialDev.fXX()+elasticStrainV/3.;
    // data.fElasticStrain[index].fYY() = factor*elasticTrialDev.fYY()+elasticStrainV/3.;
    // data.fElasticStrain[index].fZZ() = factor*elasticTrialDev.fZZ()+elasticStrainV/3.;
    // data.fElasticStrain[index].fXY() = factor*elasticTrialDev.fXY()*2.;









    //Box 7.3
    auto hydrostatic = Stress.Hydrostatic();
    //Flow vector with trial deviatory stress
    fFlowVector = FlowVector(fDeviatory);

    fDeviatory *= (1. - data.fPlasticMultiplier[index] * 3. * fShearModulus / fVonMisesStress);
    Stress = fDeviatory + hydrostatic;
    
    Tensor epsilonUpdated(fDeviatory);
    epsilonUpdated /= (2. * fShearModulus);
    epsilonUpdated.fXY() *= 2.;
    epsilonUpdated.fXZ() *= 2.;
    epsilonUpdated.fYZ() *= 2.;
    
    
    Tensor Ident;
    double epslion_e_trial = (data.fElasticStrain[index]).Trace();
    Ident.Identity();
    Ident *= epslion_e_trial / 3.;
    epsilonUpdated += Ident;
    //Not needed, already tested comparing with 3d
    // if (fElasticModel->Dimension() == 2 && !fPlaneStress){
    //     epsilonUpdated.fXZ() = 0.;
    //     epsilonUpdated.fYZ() = 0.;
    //     epsilonUpdated.fZZ() = 0.;
    // }
    //Update elastic strain
    data.fElasticStrain[index] = epsilonUpdated;

    // data.fPlasticStrain[index] += data.fPlasticMultiplier[index];
    
}   