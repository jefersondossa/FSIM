#include "DruckerPrager.h"

DruckerPrager::DruckerPrager(WeakForm *elast, REAL phi, REAL psi, bool oe) : PlasticityModel(elast){

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
    REAL tanphi = tan(fInternalFriction);
    REAL tanpsi = tan(fDilatancyAngle);
    fEta = 3. * tanphi / sqrt(9. + 12.*tanphi*tanphi);
    fXi = 3. / sqrt(9. + 12.*tanphi*tanphi);
    fEtaBar = 3. * tanpsi / sqrt(9. + 12.*tanpsi*tanpsi);


    fAlpha = fXi/fEta;
    fBeta = fXi/fEtaBar;
    //Equation 8.120

}

void DruckerPrager::ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor3D &Stress){
    
    MatrixDouble fTangentTensor(6,6);
    fTangentTensor.setZero();
    //std::cout << "Put the tangent elastoplastic tensor here.\n";

    if (fPlaneStress){
        PanicButton();
    } else {
        //Elastic Tensor
        fTangentTensor = 2.*fShearModulus*fIdentity4Dev + fBulkModulus*fId2xId2;
        // if (fApex){
        //     fTangentTensor3D = fBulkModulus*(1.-fBulkModulus/(fBulkModulus+fAlpha*fBeta*fHardening))*fId2xId2;  
        // }else{
        //     REAL A = 1./(fShearModulus + fBulkModulus*fEta*fEtaBar + fXi*fXi*fHardening);
        //     REAL sq2 = sqrt(2.);
        //     REAL devstrainnorm = fTrialDevStrain.Norm();
        //     if (fabs(devstrainnorm) > 1.e-10){
        //         fTrialDevStrain /= devstrainnorm;
        //     } else {
        //         fTrialDevStrain.Zero();
        //     }
        //     // fTrialDevStrain.Zero();
        //     Tensor3D Ident;
        //     Ident.Identity();
        //     auto DxD = fTrialDevStrain.TensorProduct(fTrialDevStrain);
        //     auto DxI = fTrialDevStrain.TensorProduct(Ident);
        //     auto IxD = Ident.TensorProduct(fTrialDevStrain);
            
        //     fTangentTensor3D = 2. * fShearModulus * (1. - data.fPlasticMultiplier[index]/(sq2 * devstrainnorm)) * fIdentity4Dev
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


REAL DruckerPrager::YieldFunction(int &index, IntPointData &data, Tensor3D &Stress){
    REAL YF = 0.;
    REAL fCohesion = 0.;
    if (fPlaneStress){
         
    } else {
        //Equation 8.101
        REAL p = Stress.Trace()/3.;
        
        fUniaxialYield(data.fEffectivePlasticStrain[index],fCohesion,fHardening);
        YF = sqrt(Stress.J2()) + fEta*p - fXi*fCohesion;
    }
    REAL tanphi = tan(fInternalFriction);
    fAlpha2 = tanphi / sqrt(9. + 12.*tanphi*tanphi);
    fK = 3.*fCohesion / sqrt(9. + 12.*tanphi*tanphi); 
    REAL YF2 = sqrt(Stress.J2()) - fAlpha2 * Stress.I1() - fK;
    return YF;
}

REAL DruckerPrager::PlasticMultiplier(int &index, IntPointData &data, Tensor3D &Stress){
    //Newton-Raphson to find plastic multiplier

    // fEtaBar criasdo no .h
    REAL plinha = sqrt(Stress.J2()) / fK * (1. + 9.*fAlpha2*fAlpha2*fBulkModulus/fShearModulus);
    REAL sigmaDepsilon = data.fElasticStrainIncrement[index].DoubleContraction(Stress);
    REAL demm = data.fElasticStrainIncrement[index].Trace();
    REAL lambda = (sigmaDepsilon-(fK/(3.*fAlpha2))*demm*(plinha-1.))/(plinha*fK);
    REAL dGamma = 0.;
    if (fPlaneStress){
        PanicButton();
    } else {
        REAL p = Stress.Trace()/3.;
        REAL sqJ2 = sqrt(Stress.J2());

        REAL fCohesion = 0.;
        fUniaxialYield(data.fEffectivePlasticStrain[index],fCohesion,fHardening);
        REAL PhiTil = sqJ2 + fEta*p -fXi*fCohesion;

        int maxiter = 10;
        int iter = 0;
        while (fabs(PhiTil) > 1.e-5){
            iter++;
            if (iter == maxiter) break;
            REAL d = -fShearModulus - fBulkModulus*fEtaBar*fEta - fHardening*fXi*fXi;
            dGamma -= PhiTil/d;

            data.fEffectivePlasticStrain[index] += fXi*dGamma;
            fUniaxialYield(data.fEffectivePlasticStrain[index],fCohesion,fHardening);
            PhiTil = sqJ2 - fShearModulus*dGamma + fEta*(p-fBulkModulus*fEtaBar*dGamma)- fXi*fCohesion;
        }

        
        if((sqJ2 - fShearModulus*dGamma)>= 0){
            fApex = false;
            return dGamma;
        } else{
            // PanicButton();
            //fAlpha and fBeta created on .h
            fApex = true;
            REAL ptrial = Stress.Trace()/3.;
            REAL fCohesion = 0.;
            fUniaxialYield(data.fEffectivePlasticStrain[index],fCohesion,fHardening);
            REAL r = fCohesion*fBeta - ptrial;
            REAL depsilon =0;
            //Box 8.10
            int maxiter = 10;
            int iter = 0;
            while (fabs(r) > 1.e-5){
                iter++;
                if (iter == maxiter) break;
                REAL d = fAlpha*fBeta*fHardening + fBulkModulus;
                depsilon -= r/d;

                data.fEffectivePlasticStrain[index] += fAlpha*depsilon;
                fUpdatedPressure = ptrial - fBulkModulus*depsilon;
                fUniaxialYield(data.fEffectivePlasticStrain[index],fCohesion,fHardening);
                r = fBeta*fCohesion - fUpdatedPressure;
            }
        }
    }
    
    
    

    return dGamma;
}

void DruckerPrager::UpdateStateVariables(int &index, IntPointData &data, Tensor3D &Stress){
    
    auto strial = Stress.Deviatory();
    Tensor3D epsilonUpdated(strial);
    Tensor3D Ident;
    
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
    fTrialDevStrain = data.fElasticStrain[index].Deviatory();
    data.fElasticStrain[index] = epsilonUpdated;

}   