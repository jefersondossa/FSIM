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


    fAlpha = fXi/fEta;
    fBeta = fXi/fEtaBar;
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
        
        fUniaxialYield(data.fPlasticStrain[index],fCohesion,fHardening);
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
    double dGamma = 0.;
    if (fPlaneStress){
        PanicButton();
    } else {
        double p = Stress.Trace()/3.;
        double sqJ2 = sqrt(Stress.J2());

        double fCohesion = 0.;
        fUniaxialYield(data.fPlasticStrain[index],fCohesion,fHardening);
        double PhiTil = sqJ2 + fEta*p -fXi*fCohesion;

        int maxiter = 10;
        int iter = 0;
        while (fabs(PhiTil) > 1.e-5){
            iter++;
            if (iter == maxiter) break;
            double d = -fShearModulus - fBulkModulus*fEtaBar*fEta - fHardening*fXi*fXi;
            dGamma -= PhiTil/d;

            data.fPlasticStrain[index] += fXi*dGamma;
            fUniaxialYield(data.fPlasticStrain[index],fCohesion,fHardening);
            PhiTil = sqJ2 - fShearModulus*dGamma + fEta*(p-fBulkModulus*fEtaBar*dGamma)- fXi*fCohesion;
        }

        
        if((sqJ2 - fShearModulus*dGamma)>= 0){
            fApex = false;
            return dGamma;
        } else{
            // PanicButton();
            //fAlpha and fBeta created on .h
            fApex = true;
            double ptrial = Stress.Trace()/3.;
            double fCohesion = 0.;
            fUniaxialYield(data.fPlasticStrain[index],fCohesion,fHardening);
            double r = fCohesion*fBeta - ptrial;
            double depsilon =0;
            //Box 8.10
            int maxiter = 10;
            int iter = 0;
            while (fabs(r) > 1.e-5){
                iter++;
                if (iter == maxiter) break;
                double d = fAlpha*fBeta*fHardening + fBulkModulus;
                depsilon -= r/d;

                data.fPlasticStrain[index] += fAlpha*depsilon;
                fUpdatedPressure = ptrial - fBulkModulus*depsilon;
                fUniaxialYield(data.fPlasticStrain[index],fCohesion,fHardening);
                r = fBeta*fCohesion - fUpdatedPressure;
            }
        }
    }


    return dGamma;
}

void DruckerPrager::UpdateStateVariables(int &index, IntPointData &data, Tensor &Stress){
    
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
    fTrialDevStrain = data.fElasticStrain[index].Deviatory();
    data.fElasticStrain[index] = epsilonUpdated;

}   