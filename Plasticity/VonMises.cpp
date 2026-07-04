#include "VonMises.h"

VonMises::VonMises(WeakForm *elast) : PlasticityModel(elast){
    fFlowVector.Zero();
    fDeviatory.Zero();

    if (fPlaneStress) {
        fMatP.resize(3,3);//Box 9.3
        fMatP.setZero();
        fMatP(0,0) = fMatP(1,1) = 2./3.;
        fMatP(0,1) = fMatP(1,0) =-1./3.;
        fMatP(2,2) = 2.;
    }
}

void VonMises::ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor3D &Stress){
    
    MatrixDouble fTangentTensor(6,6);

    fTangentTensor.setZero();
    
    MatrixDouble NxN = fFlowVector.TensorProduct(fFlowVector);
    
    REAL Afactor = 2. * fShearModulus *(1.-3.*fShearModulus*data.fPlasticMultiplier[index]/fVonMisesStress);
    REAL Bfactor = 6. * fShearModulus * fShearModulus * (data.fPlasticMultiplier[index]/fVonMisesStress - 1./(3.*fShearModulus+fHardening));

    //Elastic operator 7.107
    if (fPlaneStress){
        // PanicButton();
        fTangentTensor = 2.*fShearModulus*fIdentity4Dev + fBulkModulus*fId2xId2;

    } else {
        //Elastic Tensor
        fTangentTensor = 2.*fShearModulus*fIdentity4Dev + fBulkModulus*fId2xId2;
        //Elastoplastic tensor
        //Continuum operator
        // fTangentTensor3D = 2.*fShearModulus*fIdentity4Dev + fBulkModulus*fId2xId2 - 6.*fShearModulus*fShearModulus/(3*fShearModulus+fHardening)*NxN;
        //Consistent operator
        // fTangentTensor3D = Afactor*fIdentity4Dev + Bfactor*NxN + fBulkModulus*fId2xId2;  
    }
    
    
    // auto elast = fElasticModel->ConstitutiveMatrix();
    // std::cout << "Elastic - \n" << elast << std::endl;
    // std::cout << "Tangent - \n" << fTangentTensor3D << std::endl;

    if (!fPlaneStress){
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
    }
    
    fElasticModel->ComputeStiffness(index,data,Stiffness);
};
    
void VonMises::ComputeError(IntPointData &data, VecDouble &errors){
    fElasticModel->ComputeError(data,errors);
};


REAL VonMises::YieldFunction(int &index, IntPointData &data, Tensor3D &Stress){
    REAL YF = 0.;
    if (fPlaneStress){
        //Box 9.3
        VecDouble Sigma(3);
        Sigma[0] = Stress.fData[0];
        Sigma[1] = Stress.fData[1];
        Sigma[2] = Stress.fData[3];
        REAL sigmay;
        fUniaxialYield(data.fEffectivePlasticStrain[index],sigmay,fHardening);
        YF = 0.5 * Sigma.transpose()*fMatP*Sigma-sigmay*sigmay/3.; 
    } else {
        //Eq (7.74)
        REAL sigmay;
        fUniaxialYield(data.fEffectivePlasticStrain[index],sigmay,fHardening);
        YF = sqrt(3.*Stress.J2()) - sigmay;
    }
    return YF;
}

Tensor3D VonMises::FlowVector(Tensor3D &Stress){
    //Eq (7.117)
    REAL devnorm = Stress.DeviatoryNorm();
    auto temp = Stress.Deviatory();
    temp *= 1. / devnorm;
    auto norm = temp.Norm();
    return temp;
}

REAL VonMises::PlasticMultiplier(int &index, IntPointData &data, Tensor3D &Stress){
    //Newton-Raphson to find plastic multiplier
    REAL DGAMA = 0.;
    if (fPlaneStress){//Box 9.5
        REAL A1 = (Stress.fXX()+Stress.fYY())*(Stress.fXX()+Stress.fYY());
        REAL A2 = (Stress.fYY()-Stress.fXX())*(Stress.fYY()-Stress.fXX());
        REAL A3 = Stress.fXY() * Stress.fXY();
        REAL XI = A1 / 6. + 0.5 * A2 + 2. * A3;
        REAL SIGMAY;
        fUniaxialYield(data.fEffectivePlasticStrain[index],SIGMAY,fHardening);
        REAL PHI = 0.5 * XI - SIGMAY*SIGMAY / 3.;
        REAL EPBARN = data.fEffectivePlasticStrain[index];
        REAL SQR2D3 = sqrt(2./3.);
        REAL isValidStress = false;
        if (PHI/SIGMAY > 1.e-6){
            REAL SQRTXI = sqrt(XI);
            REAL B1 = 1.;
            REAL B2 = 1.;
            REAL FMODU = fYoungModulus/(3.*(1.-fPoissonRatio));
            int NITER = 0;
            REAL EPBAR = EPBARN;
            
            while (NITER < 100){
                NITER ++;
                //Compute residual derivative
                REAL HSLOPE = fHardening;
                REAL DXI = -A1*FMODU/(3.*B1*B1*B1)-2.*fShearModulus*(A2+4.*A3)/(B2*B2*B2);
                REAL HBAR = 2.*SIGMAY*HSLOPE*SQR2D3*(SQRTXI+DGAMA*DXI/(2.*SQRTXI));
                REAL DPHI=0.5*DXI-HBAR/3.;
                //Compute Newton-Raphson increment and update equation variable DGAMA
                DGAMA=DGAMA-PHI/DPHI;
                //Compute new residual (yield function value)
                B1=1.+FMODU*DGAMA;
                B2=2.+2.*fShearModulus*DGAMA;
                XI=(1./6.)*A1/(B1*B1)+(0.5*A2+2.*A3)/(B2*B2);
                SQRTXI=sqrt(XI);
                EPBAR=EPBARN+DGAMA*SQR2D3*SQRTXI;
                fUniaxialYield(EPBAR,SIGMAY,fHardening);
                PHI=0.5*XI-SIGMAY*SIGMAY/3.;
                //Check for convergence
                REAL RESNOR=fabs(PHI/SIGMAY);
                if (RESNOR <= 1.e-6){
                    isValidStress = true;
                    // update accumulated plastic strain
                    data.fEffectivePlasticStrain[index] = EPBAR;
                    REAL ASTAR1=3.*(1.-fPoissonRatio)/(3.*(1.-fPoissonRatio)+fYoungModulus*DGAMA);
                    REAL ASTAR2=1./(1.+2.*fShearModulus*DGAMA);
                    REAL A11=0.5*(ASTAR1+ASTAR2);
                    REAL A22 = A11;
                    REAL A12=0.5*(ASTAR1-ASTAR2);
                    REAL A21 = A12;
                    REAL A33 = ASTAR2;
                    VecDouble StressUp(3);
                    StressUp[0] = A11*Stress.fXX()+A12*Stress.fYY();
                    StressUp[1] = A21*Stress.fXX()+A22*Stress.fYY();
                    StressUp[2] = A33*Stress.fXY();
                    Stress.fXX() = StressUp[0];
                    Stress.fYY() = StressUp[1];
                    Stress.fXY() = StressUp[2];
                    // compute corresponding elastic (engineering) strain components
                    REAL FACTG=1./(2.*fShearModulus);
                    REAL P=(StressUp[0]+StressUp[1])/3.;
                    REAL EEV=P/fBulkModulus;
                    REAL EEVD3=EEV/3.;
                    VecDouble ElasticStrainUp(4);
                    ElasticStrainUp[0] = FACTG*(2.*StressUp[0]/3.-StressUp[1]/3.)+EEVD3;
                    ElasticStrainUp[1] = FACTG*(2.*StressUp[1]/3.-StressUp[0]/3.)+EEVD3;
                    ElasticStrainUp[2] = FACTG*StressUp[2]*2.;
                    ElasticStrainUp[3] = -fPoissonRatio/(1.-fPoissonRatio)*(ElasticStrainUp[0]+ElasticStrainUp[1]);
                    data.fElasticStrain[index].fXX() = ElasticStrainUp[0];
                    data.fElasticStrain[index].fYY() = ElasticStrainUp[1];
                    data.fElasticStrain[index].fXY() = ElasticStrainUp[2];
                    data.fElasticStrain[index].fZZ() = ElasticStrainUp[3];
                    break;
                }
            }

        }
        if (!isValidStress){
            PanicButton();
        }
        


        // REAL splus2 = (Stress.fXX()+Stress.fYY())*(Stress.fXX()+Stress.fYY());
        // REAL smin2 = (Stress.fYY()-Stress.fXX())*(Stress.fYY()-Stress.fXX());
        // REAL t2 = Stress.fXY() * Stress.fXY();
        // REAL xi = splus2 / 6. + smin2 / 2 + 2. * t2;
        // REAL sigmay;
        // fUniaxialYield(data.fEffectivePlasticStrain[index],sigmay,fHardening);
        // REAL PhiTil = 0.5 * xi - sigmay*sigmay / 3.;
        // while (fabs(PhiTil) > 1.e-5){
        //     REAL dXi = (- splus2 * fYoungModulus/(1.-fPoissonRatio)) / (9.*pow(1.+fYoungModulus*dGamma/(3.*(1.-fPoissonRatio)),3))
        //                - 2.*fShearModulus * (smin2+4.*t2) / pow(1.+2.*fShearModulus*dGamma,3);
        //     fUniaxialYield(data.fEffectivePlasticStrain[index]+dGamma*sqrt(2.*xi/3.),sigmay,fHardening);
        //     REAL fHardBar = 2.*sigmay*fHardening*sqrt(2./3.)*(sqrt(xi)+dGamma*dXi/(2.*sqrt(xi)));
        //     REAL dPhiTil = 0.5*dXi-fHardBar/3.;
        //     dGamma -= PhiTil/dPhiTil;
        //     xi = splus2/(6.*pow(1.+fYoungModulus*dGamma/(3.*(1.-fPoissonRatio)),2))+
        //         (0.5*smin2+2.*t2)/pow(1.+2.*fShearModulus*dGamma,2);
        //     fUniaxialYield(data.fEffectivePlasticStrain[index]+dGamma*sqrt(2.*xi/3.),sigmay,fHardening);
        //     PhiTil = 0.5 * xi - sigmay*sigmay / 3.;
        // }
        // data.fEffectivePlasticStrain[index] += dGamma * sqrt(2.*xi/3.);
    } else { //Box 7.4 Souza Neto
        fDeviatory = Stress.Deviatory();
        fVonMisesStress = sqrt(1.5*fDeviatory.DoubleContraction(fDeviatory))+3.*data.fPlasticMultiplier[index];
        REAL sigmay = 0.;
        fUniaxialYield(data.fEffectivePlasticStrain[index],sigmay,fHardening);
        REAL PhiTil = fVonMisesStress-sigmay;
        
        while (fabs(PhiTil/sigmay) > 1.e-5){
            REAL d = -3. * fShearModulus - fHardening;
            DGAMA -= PhiTil/d;
            fUniaxialYield(data.fEffectivePlasticStrain[index]+DGAMA,sigmay,fHardening);
            PhiTil = fVonMisesStress - 3.*fShearModulus*DGAMA-sigmay;
        }
        data.fEffectivePlasticStrain[index] += DGAMA;
    }

    return DGAMA;
}

void VonMises::UpdateStateVariables(int &index, IntPointData &data, Tensor3D &Stress){
    
    if (fPlaneStress){
        return;
        // MatrixDouble MatA(3,3);
        // MatA.setZero();
        // REAL A11 = 3.*(1.-fPoissonRatio)/(3.*(1.-fPoissonRatio)+fYoungModulus*data.fPlasticMultiplier[index]);
        // REAL A22 = 1./(1.+2.*fShearModulus*data.fPlasticMultiplier[index]);
        // REAL A33 = A22;
        // MatA(0,0) = MatA(1,1) = 0.5 * (A11+A22);
        // MatA(1,0) = MatA(0,1) = 0.5 * (A11-A22);
        // MatA(2,2) = A33;
        // VecDouble auxStress(3);
        // auxStress[0] = Stress.fXX();
        // auxStress[1] = Stress.fYY();
        // auxStress[2] = Stress.fXY();
        // VecDouble res = MatA * auxStress;

        // Stress.Zero();
        // Stress.fXX() = res[0];
        // Stress.fYY() = res[1];
        // Stress.fXY() = res[2];
        // //Eq. 9.4
        // REAL alpha = (3.*fBulkModulus - 2.*fShearModulus) / (3.*fBulkModulus + 4.*fShearModulus);
        // MatrixDouble MatD(3,3);
        // MatD.setZero();
        // MatD(0,0) = MatD(1,1) = 1. + alpha;
        // MatD(0,1) = MatD(1,0) = alpha;
        // MatD(2,2) = 0.5;
        // MatD *= 2.*fShearModulus;
        
        // data.fElasticStrain[index].Zero();
        // VecDouble elastStrain = MatD.inverse()*res;
        // data.fElasticStrain[index].fXX() = elastStrain[0];
        // data.fElasticStrain[index].fYY() = elastStrain[1];
        // data.fElasticStrain[index].fXY() = elastStrain[2];
        // data.fElasticStrain[index].fZZ() = -fPoissonRatio/(1.-fPoissonRatio) * (elastStrain[0]+elastStrain[1]);

    } else {
        //Box 7.3
        auto hydrostatic = Stress.Hydrostatic();
        //Flow vector with trial deviatory stress
        fFlowVector = FlowVector(fDeviatory);

        fDeviatory *= (1. - data.fPlasticMultiplier[index] * 3. * fShearModulus / fVonMisesStress);
        Stress = fDeviatory + hydrostatic;
        
        Tensor3D epsilonUpdated(fDeviatory);
        epsilonUpdated /= (2. * fShearModulus);
        epsilonUpdated.fXY() *= 2.;
        epsilonUpdated.fXZ() *= 2.;
        epsilonUpdated.fYZ() *= 2.;
        
        
        Tensor3D Ident;
        REAL epslion_e_trial = (data.fElasticStrain[index]).Trace();
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

        // data.fEffectivePlasticStrain[index] += data.fPlasticMultiplier[index];
    }
}   