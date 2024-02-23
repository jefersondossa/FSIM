#include "Tresca.h"
#include "ElasticTruss.h"
#include "Elasticity2D.h"
#include "PositionalTruss.h"
#include "ElasticityPositional2D.h"

Tresca::Tresca(WeakForm *elast) : PlasticityModel(elast){

   
}


void Tresca::ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor &Stress){
    
    MatrixDouble fTangentTensor(6,6);
    fTangentTensor.setZero();
    MatrixDouble dSdEpsilon(3,3);
    dSdEpsilon.setZero();

    switch (fReturnDirection)
    {
    case EMainPlane:
        {
            //Eq. 8.41
            double f = 2.*fShearModulus / (4.*fShearModulus + fHardening);
            dSdEpsilon(0,0) = dSdEpsilon(2,2) = 2.*fShearModulus * (1. - f);
            dSdEpsilon(1,1) = 2.*fShearModulus;
            dSdEpsilon(0,2) = dSdEpsilon(2,0) = 2.*fShearModulus * f;
        }
        break;
    case ERightCorner:
        {
            //Eq. 8.52
            MatrixDouble d(2,2);
            d(0,0) = d(1,1) = -4.*fShearModulus-fHardening;
            d(0,1) = d(1,0) = -2.*fShearModulus-fHardening;
            double detD = d.determinant();

            dSdEpsilon(0,0) = 2.*fShearModulus * (1. - 8.*fShearModulus*fShearModulus/detD);
            dSdEpsilon(1,1) = 2.*fShearModulus * (1. + 2.*fShearModulus*d(0,0)/detD);
            dSdEpsilon(2,2) = 2.*fShearModulus * (1. + 2.*fShearModulus*d(1,1)/detD);
            dSdEpsilon(1,0) = dSdEpsilon(2,0) = 8.*fShearModulus*fShearModulus*fShearModulus/detD;
            dSdEpsilon(0,1) = 4.*fShearModulus*fShearModulus*(d(0,1)-d(0,0))/detD;
            dSdEpsilon(0,2) = 4.*fShearModulus*fShearModulus*(d(1,0)-d(1,1))/detD;
            dSdEpsilon(1,2) = -4.*fShearModulus*fShearModulus*(d(1,0))/detD;
            dSdEpsilon(2,1) = -4.*fShearModulus*fShearModulus*(d(0,1))/detD;
        }
        break;
    case ELeftCorner:
        {
            //Eq. 8.53
            MatrixDouble d(2,2);
            d(0,0) = d(1,1) = -4.*fShearModulus-fHardening;
            d(0,1) = d(1,0) = -2.*fShearModulus-fHardening;
            double detD = d.determinant();

            dSdEpsilon(2,2) = 2.*fShearModulus * (1. - 8.*fShearModulus*fShearModulus/detD);
            dSdEpsilon(1,1) = 2.*fShearModulus * (1. + 2.*fShearModulus*d(0,0)/detD);
            dSdEpsilon(0,0) = 2.*fShearModulus * (1. + 2.*fShearModulus*d(1,1)/detD);
            dSdEpsilon(0,2) = dSdEpsilon(1,2) = 8.*fShearModulus*fShearModulus*fShearModulus/detD;
            dSdEpsilon(2,1) = 4.*fShearModulus*fShearModulus*(d(0,1)-d(0,0))/detD;
            dSdEpsilon(2,0) = 4.*fShearModulus*fShearModulus*(d(1,0)-d(1,1))/detD;
            dSdEpsilon(1,0) = -4.*fShearModulus*fShearModulus*(d(1,0))/detD;
            dSdEpsilon(0,1) = -4.*fShearModulus*fShearModulus*(d(0,1))/detD;
        }
        break;
    default:
        PanicButton();
        break;
    }

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                double dkj = k == j ? 1. : 0;
                fTangentTensor(i,j) = dSdEpsilon(i,k)*(dkj - 1./3.) + fBulkModulus;
            }
        }
    }
    

    // std::cout << "Put the tangent elastoplastic tensor here.\n";
    // PanicButton();   

    if (fElasticModel->Dimension() == 2){
        MatrixDouble fTangent2D(3,3);
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                fTangent2D(i,j) = fTangentTensor(i,j);
            }
        }
        fElasticModel->ConstitutiveMatrix() = fTangent2D;
    } else if (fElasticModel->Dimension() == 3) {
        std::cout << "Please check this implementation" << std::endl;
        PanicButton();
        fElasticModel->ConstitutiveMatrix() = fTangentTensor;
    } else {
        PanicButton();
    }
    
    fElasticModel->ComputeStiffness(index,data,Stiffness);
};

double Tresca::YieldFunction(int &index, IntPointData &data, Tensor &Stress){
    double YF = 0.;

    if (fPlaneStress){
         
    } else {
        VecDouble eigenvalues;
        auto deviatory = Stress.Deviatory();
        deviatory.SpectralDecomposition(eigenvalues,fPrincipalDirections);
        fS1 = eigenvalues[0];
        fS2 = eigenvalues[1];
        fS3 = eigenvalues[2];
        double sigmay = 0.;
        fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
        YF = fS1-fS3-sigmay;
    }
    return YF;
}

double Tresca::PlasticMultiplier(int &index, IntPointData &data, Tensor &Stress){
    //Newton-Raphson to find plastic multiplier
    double dGamma = 0.;
    double s1 = 0.;
    double s2 = 0.;
    double s3 = 0.;

    if (fPlaneStress){
        
    } else {
        // First attempt - return to main plane
        double sigmay = 0.;
        fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
        double PhiTil = fS1-fS3-sigmay;
        while (fabs(PhiTil) > 1.e-8){
            double d = -4. * fShearModulus - fHardening;
            dGamma -= PhiTil/d;
            fUniaxialYield(data.fPlasticStrain[index]+dGamma,sigmay,fHardening);
            PhiTil = fS1 - fS3 - 4.*fShearModulus*dGamma - sigmay;
        }
        s1 = fS1 - 2. * fShearModulus * dGamma;
        s2 = fS2;
        s3 = fS3 + 2. * fShearModulus * dGamma;
        fReturnDirection = EMainPlane;
        fS1 = s1;
    fS3 = s2;
    fS3 = s3;

    return dGamma;

        //Check validity of main plane return
        if (s1 >= s2 && s2 >= s3){
            //Is in the main plane
            data.fPlasticStrain[index] += dGamma;
        } else {
            std::cout << "Probably there is a bug here. Please debug me." << std::endl;
            PanicButton();
            //Return to corner - Box 8.3
            VecDouble PhiTilAB(2);
            double sa = fS1-fS3;
            double sb;
            fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
            if (fS1 + fS3 - 2.*fS2 > 0){//Apply return to the RIGHT corner
                fReturnDirection = ERightCorner;
                sb = fS1-fS2;
            } else {//Apply return to the LEFT corner
                fReturnDirection = ELeftCorner;
                sb = fS2-fS3;
            }
            PhiTilAB[0] = sa-sigmay;
            PhiTilAB[1] = sb-sigmay;

            VecDouble dGammaAB(2);
            dGammaAB.setZero();
            MatrixDouble d(2,2);
            d(0,0) = d(1,1) = -4.*fShearModulus - fHardening;
            d(0,1) = d(1,0) = -2.*fShearModulus - fHardening;
            MatrixDouble dinv = d.inverse();
            int niterations = 0;
            double epn = data.fPlasticStrain[index];
            while (fabs(PhiTilAB[0])+fabs(PhiTilAB[1]) > 1.e-8){
                niterations++;
                double dGammaBarra = dGammaAB[0] + dGammaAB[1];
                data.fPlasticStrain[index] += dGammaBarra;
                dGammaAB -= dinv * PhiTilAB;
                fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
                PhiTilAB[0] = sa - 2.*fShearModulus*(2.*dGammaAB[0]+dGammaAB[1])-sigmay;
                PhiTilAB[1] = sb - 2.*fShearModulus*(dGammaAB[0]+2.*dGammaAB[1])-sigmay;
            }
            
            dGamma=dGammaAB[0]+dGammaAB[1];
             
            if (fS1 + fS3 - 2.*fS2 > 0){//Apply return to the RIGHT corner
                s1 = fS1 - 2. * fShearModulus * (dGammaAB[0]+dGammaAB[1]);
                s2 = fS2 + 2. * fShearModulus * dGammaAB[1];
                s3 = fS3 + 2. * fShearModulus * dGammaAB[0];
            } else {//Apply return to the LEFT corner
                s1 = fS1 - 2. * fShearModulus * dGammaAB[0];
                s2 = fS2 + 2. * fShearModulus * dGammaAB[1];
                s3 = fS3 + 2. * fShearModulus * (dGammaAB[0]+dGammaAB[1]);
            }  
        }
    }
    fS1 = s1;
    fS3 = s2;
    fS3 = s3;

    return dGamma;
}

void Tresca::UpdateStateVariables(int &index, IntPointData &data, Tensor &Stress){
    //Box 8.1
    double pressure = Stress.Trace() / 3.;

    MatrixDouble sn1(3,3),aux1,aux2,aux3;
    // aux1= fPrincipalDirections.col(0) * fPrincipalDirections.col(0).transpose();
    // aux2= fPrincipalDirections.col(1) * fPrincipalDirections.col(1).transpose();
    // aux3= fPrincipalDirections.col(2) * fPrincipalDirections.col(2).transpose();
    // std::cout << "S1 = \n" << aux1 << std::endl;
    // std::cout << "S2 = \n" << aux2 << std::endl;
    // std::cout << "S3 = \n" << aux3 << std::endl;
    sn1 = (pressure + fS1) * fPrincipalDirections.col(0) * fPrincipalDirections.col(0).transpose()
         +(pressure + fS2) * fPrincipalDirections.col(1) * fPrincipalDirections.col(1).transpose()
         +(pressure + fS3) * fPrincipalDirections.col(2) * fPrincipalDirections.col(2).transpose();
    Tensor StressN1(sn1);
    Stress = StressN1;

    Tensor epsilonUpdated(StressN1.Deviatory());
    epsilonUpdated /= (2. * fShearModulus);
    Tensor Ident;
    double epslion_e_trial = data.fElasticStrain[index].Trace() / 3.;
    Ident.Identity();
    Ident *= epslion_e_trial;
    epsilonUpdated += Ident;

}   