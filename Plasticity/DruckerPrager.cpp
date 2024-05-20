#include "DruckerPrager.h"

DruckerPrager::DruckerPrager(WeakForm *elast, double phi, double psi, bool oe) : PlasticityModel(elast){

    fInternalFriction = phi;
    fOuterEdges = oe;
    fDilatancyAngle = psi;

    if (fOuterEdges){
        //Equation 6.122 and 6.123
        fEta = 6.*sin(fInternalFriction)/(sqrt(3.)*(3.-sin(fInternalFriction)));
        fXi = 6.*cos(fInternalFriction)/(sqrt(3.)*(3.-sin(fInternalFriction)));
        fEtaBar = 6.*sin(fDilatancyAngle)/(sqrt(3.)*(3.-sin(fDilatancyAngle)));
        //fEtaBar = fEta com traço - Equation 6.163
    }
    else{
        fEta = 6.*sin(fInternalFriction)/(sqrt(3.)*(3.+sin(fInternalFriction)));
        fXi = 6.*cos(fInternalFriction)/(sqrt(3.)*(3.+sin(fInternalFriction))); 
        fEtaBar = 6.*sin(fDilatancyAngle)/(sqrt(3.)*(3.+sin(fDilatancyAngle)));
    }
    double tanphi = tan(fInternalFriction);
    double tanpsi = tan(fDilatancyAngle);
    fEta = 3. * tanphi / sqrt(9. + 12.*tanphi*tanphi);
    fXi = 3. / sqrt(9. + 12.*tanphi*tanphi);
    fEtaBar = 3. * tanpsi / sqrt(9. + 12.*tanpsi*tanpsi);


    fAlpha = fXi/fEtaBar;
    fBeta = fXi/fEta;
    //Equation 8.120

}

void DruckerPrager::ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor &Stress){
    
    MatrixDouble fTangentTensor(6,6);
    fTangentTensor.setZero();
    //std::cout << "Put the tangent elastoplastic tensor here.\n";

    if (fPlaneStress){
        PanicButton();
    } else {
        //Elastic Tensor
        fTangentTensor = 2.*fShearModulus*fIdentity4Dev + fBulkModulus*fId2xId2;
        // if (fApex){
        //     fTangentTensor = fBulkModulus*(1.-fBulkModulus/(fBulkModulus+fAlpha*fBeta*fHardening))*fId2xId2;  
        // }else{
        //     double A = 1./(fShearModulus + fBulkModulus*fEta*fEtaBar + fXi*fXi*fHardening);
        //     double sq2 = sqrt(2.);
        //     double devstrainnorm = fTrialDevStrain.Norm();
        //     if (fabs(devstrainnorm) > 1.e-10){
        //         fTrialDevStrain /= devstrainnorm;
        //     } else {
        //         fTrialDevStrain.Zero();
        //     }
        //     // fTrialDevStrain.Zero();
        //     Tensor Ident;
        //     Ident.Identity();
        //     auto DxD = fTrialDevStrain.TensorProduct(fTrialDevStrain);
        //     auto DxI = fTrialDevStrain.TensorProduct(Ident);
        //     auto IxD = Ident.TensorProduct(fTrialDevStrain);
            
        //     fTangentTensor = 2. * fShearModulus * (1. - data.fPlasticMultiplier[index]/(sq2 * devstrainnorm)) * fIdentity4Dev
        //                    + 2. * fShearModulus * (data.fPlasticMultiplier[index]/(sq2 * devstrainnorm) - fShearModulus*A) * DxD
        //                    - sq2*fShearModulus*A*fBulkModulus*(fEta*DxI + fEtaBar*IxD)
        //                    + fBulkModulus * (1. - fBulkModulus * fEta * fEtaBar * A)*fId2xId2;
        // }
    }

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
    
void DruckerPrager::ComputeError(IntPointData &data, VecDouble &errors){
    fElasticModel->ComputeError(data,errors);
};


double DruckerPrager::YieldFunction(int &index, IntPointData &data, Tensor &Stress){
    double YF = 0.;
    double fCohesion = 0.;
    if (fPlaneStress){
         
    } else {
        //Equation 8.101
        double p = Stress.Trace()/3.;
        
        fUniaxialYield(data.fEffectivePlasticStrain[index],fCohesion,fHardening);
        YF = sqrt(Stress.J2()) + fEta*p - fXi*fCohesion;
    }
    // double tanphi = tan(fInternalFriction);
    // fAlpha2 = tanphi / sqrt(9. + 12.*tanphi*tanphi);
    // fK = 3.*fCohesion / sqrt(9. + 12.*tanphi*tanphi); 
    // double YF2 = sqrt(Stress.J2()) - fAlpha2 * Stress.I1() - fK;
    return YF;
}

double DruckerPrager::PlasticMultiplier(int &index, IntPointData &data, Tensor &Stress){
    //Newton-Raphson to find plastic multiplier

    // fEtaBar criasdo no .h
    // double plinha = sqrt(Stress.J2()) / fK * (1. + 9.*fAlpha2*fAlpha2*fBulkModulus/fShearModulus);
    // double sigmaDepsilon = data.fElasticStrainIncrement[index].DoubleContraction(Stress);
    // double demm = data.fElasticStrainIncrement[index].Trace();
    // double lambda = (sigmaDepsilon-(fK/(3.*fAlpha2))*demm*(plinha-1.))/(plinha*fK);
    //  auto plaststrain = data.fEffectivePlasticStrain[index];
    // double dGamma = 0.;
    // if (fPlaneStress){
    //     PanicButton();
    // } else {
    //     double p = Stress.Trace()/3.;
    //     double sqJ2 = sqrt(Stress.J2());

    //     double fCohesion = 0.;
    //     fUniaxialYield(data.fEffectivePlasticStrain[index],fCohesion,fHardening);
    //     double PhiTil = sqJ2 + fEta*p -fXi*fCohesion;

    //     int maxiter = 100;
    //     int iter = 0;
    //     while (fabs(PhiTil) > 1.e-5){
    //         iter++;
    //         if (iter == maxiter) break;
    //         double d = -fShearModulus - fBulkModulus*fEtaBar*fEta - fHardening*fXi*fXi;
    //         dGamma -= PhiTil/d;

    //         data.fEffectivePlasticStrain[index] = plaststrain + fXi*dGamma;
    //         fUniaxialYield(data.fEffectivePlasticStrain[index],fCohesion,fHardening);
    //         PhiTil = sqJ2 - fShearModulus*dGamma + fEta*(p-fBulkModulus*fEtaBar*dGamma)- fXi*fCohesion;
    //     }
    //     if (iter == maxiter){
    //         std::cout << "Problem converging surface\n";
    //         // PanicButton();
    //     }
        
    //     if((sqJ2 - fShearModulus*dGamma)>= 0){
    //         fApex = false;
    //         return dGamma;
    //     } else{
    //         // PanicButton();
    //         //fAlpha and fBeta created on .h
    //         fApex = true;
    //         double ptrial = Stress.Trace()/3.;
    //         double fCohesion = 0.;
    //         fUniaxialYield(data.fEffectivePlasticStrain[index],fCohesion,fHardening);
    //         double r = fCohesion*fBeta - ptrial;
    //         double depsilon =0;
    //         //Box 8.10
    //         int maxiter = 100;
    //         int iter = 0;
           
    //         while (fabs(r) > 1.e-5){
    //             iter++;
    //             if (iter == maxiter) break;
    //             double d = fAlpha*fBeta*fHardening + fBulkModulus;
    //             depsilon -= r/d;

    //             data.fEffectivePlasticStrain[index] = plaststrain + fAlpha*depsilon;
    //             fUpdatedPressure = ptrial - fBulkModulus*depsilon;
    //             fUniaxialYield(data.fEffectivePlasticStrain[index],fCohesion,fHardening);
    //             r = fBeta*fCohesion - fUpdatedPressure;
    //         }
    //         if (iter == maxiter){
    //             std::cout << "Problem converging apex\n";
    //             PanicButton();
    //         }
    //     }
    // }

    double dGamma = 0.;
    double EETV = data.fElasticStrain[index].Trace();
    double EPBARN = data.fEffectivePlasticStrain[index];
    double PT = EETV * fBulkModulus;
    double EEVD3=EETV / 3.;
    Tensor strial;
    strial.Zero();
    double P,EPBAR;

    strial.fXX() = 2. * fShearModulus * (data.fElasticStrain[index].fXX() - EEVD3);
    strial.fYY() = 2. * fShearModulus * (data.fElasticStrain[index].fYY() - EEVD3);
    strial.fZZ() = 2. * fShearModulus * (data.fElasticStrain[index].fZZ() - EEVD3);
    strial.fXY() = 2. * fShearModulus * (data.fElasticStrain[index].fXY()*0.5);

    double VARJ2T = strial.fXY()*strial.fXY() + .5*(strial.fXX()*strial.fXX()+strial.fYY()*strial.fYY()+strial.fZZ()*strial.fZZ());
    
    double fCohesion = 0.;
    fUniaxialYield(data.fEffectivePlasticStrain[index],fCohesion,fHardening);
        
    double SQRJ2T=sqrt(VARJ2T);
    double PHI=SQRJ2T+fEta*PT-fXi*fCohesion;
    double RES=PHI;
    if (fCohesion != 0.)RES=RES/fabs(fCohesion);

    double tol = 1.e-8;
    int maxiter = 100;
    double FACTOR = 0.;

    if (RES > tol){
        fApex = false;
        int iter = 0;
        while (iter < maxiter){
            double DENOM = -fShearModulus-fBulkModulus*fEtaBar*fEta-fXi*fXi*fHardening;
            double DDGAMA = -PHI/DENOM;
            double dGamma = dGamma + DDGAMA;
            EPBAR = EPBARN + fXi*dGamma;
            fUniaxialYield(EPBAR,fCohesion,fHardening);
            double SQRJ2 = SQRJ2T - fShearModulus*dGamma;
            P = PT-fBulkModulus*fEtaBar*dGamma;
            PHI = SQRJ2+fEta*P-fXi*fCohesion;

            double RESNOR = fabs(PHI);
            
            if (fCohesion != 0) RESNOR = RESNOR/fabs(fCohesion);
            if (RESNOR <= tol){
                if (SQRJ2 == 0){
                    FACTOR = 0.;
                }else{
                    FACTOR = 1.-fShearModulus*dGamma/SQRJ2T;
                    break;
                }
            } else {
                fApex = true;
                if (fEta == 0 || fEtaBar == 0){
                    PanicButton();
                }
                double DEPV = 0.;
                EPBAR = EPBARN;
                fUniaxialYield(EPBAR,fCohesion,fHardening);
                RES = fBeta * fCohesion - PT;
                while (iter < maxiter){
                    DENOM=fAlpha*fBeta*fHardening+fBulkModulus;
                    double DDEPV=-RES/DENOM;
                    DEPV=DEPV+DDEPV;
                    EPBAR=EPBARN+fAlpha*DEPV;
                    fUniaxialYield(EPBAR,fCohesion,fHardening);
                    P = PT-fBulkModulus*DEPV;
                    RES = fBeta*fCohesion-P;
                    RESNOR = fabs(RES);
                    if (fCohesion != 0) RESNOR = RESNOR/fabs(fCohesion);
                    if (RESNOR <= tol){
                        dGamma = DEPV/fEtaBar;
                        FACTOR = 0.;
                    }

                }
            }
        }
    }

    Stress.fXX() = FACTOR * strial.fXX() + P;
    Stress.fYY() = FACTOR * strial.fYY() + P;
    Stress.fZZ() = FACTOR * strial.fZZ() + P;
    Stress.fXY() = FACTOR * strial.fXY();
    
    data.fEffectivePlasticStrain[index] = EPBAR;

    FACTOR = FACTOR/(2.*fShearModulus);
    EEVD3 = P/(fBulkModulus*3.);

    data.fElasticStrain[index].fXX() = FACTOR * strial.fXX() + EEVD3;
    data.fElasticStrain[index].fYY() = FACTOR * strial.fYY() + EEVD3;
    data.fElasticStrain[index].fZZ() = FACTOR * strial.fZZ() + EEVD3;
    data.fElasticStrain[index].fXY() = FACTOR * strial.fXY()*2.;

    return dGamma;
}

void DruckerPrager::UpdateStateVariables(int &index, IntPointData &data, Tensor &Stress){
    
    return;
    auto strial = Stress.Deviatory();
    Tensor epsilonUpdated(strial);
    Tensor Ident;
    
    if(fApex){
        Ident.Identity();
        Ident *= fUpdatedPressure;
        Stress = Ident;
    } else {
        strial *= (1. - fShearModulus*data.fPlasticMultiplier[index]/sqrt(Stress.J2()));
        fUpdatedPressure = Stress.Trace()/3. - fBulkModulus*fEtaBar*data.fPlasticMultiplier[index];
        
        Ident.Identity();
        Ident*=fUpdatedPressure;
        Stress = strial + Ident;
    }
    
    epsilonUpdated *= (1. - fShearModulus*data.fPlasticMultiplier[index]/sqrt(Stress.J2()))/(2*fShearModulus);
    epsilonUpdated.fXY() *= 2.;
    epsilonUpdated.fXZ() *= 2.;
    epsilonUpdated.fYZ() *= 2.;
    Ident.Identity();
    Ident *= fUpdatedPressure/(3.*fBulkModulus);
    epsilonUpdated += Ident;
    // epsilonUpdated.fXZ() = 0.;
    // epsilonUpdated.fYZ() = 0.;
    // epsilonUpdated.fZZ() = 0.;
    // fTrialDevStrain = data.fElasticStrain[index].Deviatory();
    data.fElasticStrain[index] = epsilonUpdated;

}   