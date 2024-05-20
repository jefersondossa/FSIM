#include "MohrCoulomb.h"
#include "ElasticTruss.h"
#include "Elasticity2D.h"

MohrCoulomb::MohrCoulomb(WeakForm *elast, double intfriction) : PlasticityModel(elast){
    fInternalFriction = intfriction;

    
}


void MohrCoulomb::ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor &Stress){
    
    MatrixDouble fTangentTensor(6,6);
    fTangentTensor.setZero();

    if (fPlaneStress){
        PanicButton();
    } else {
        //Elastic Tensor
        fTangentTensor = 2.*fShearModulus*fIdentity4Dev + fBulkModulus*fId2xId2;
        
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
    
void MohrCoulomb::ComputeError(IntPointData &data, VecDouble &errors){
    fElasticModel->ComputeError(data,errors);
};


double MohrCoulomb::YieldFunction(int &index, IntPointData &data, Tensor &Stress){
    double YF = 0.;
    double j2 = Stress.J2();
    double j3 = Stress.J3();
    double theta = -1./3. * asin(1.5*sqrt(3.)*j3/sqrt(j2*j2*j2));
#ifdef DEBUG_BUILD
    if (theta < -M_PI/6. || theta > M_PI/6.){
        PanicButton();
    }
#endif
    double i1 = Stress.I1();
    double fCohesion = 0.;
    fUniaxialYield(data.fEffectivePlasticStrain[index],fCohesion,fHardening);

    YF = i1*sin(fInternalFriction) + sqrt(j2)*cos(theta) 
       - sqrt(j2) * sin(theta) * sin(fInternalFriction) / 3.
       - fCohesion * cos(fInternalFriction);

    return YF;
}

double MohrCoulomb::PlasticMultiplier(int &index, IntPointData &data, Tensor &Stress){
    //Newton-Raphson to find plastic multiplier

    // fEtaBar criasdo no .h
    double dGamma;
    return dGamma;
}

void MohrCoulomb::UpdateStateVariables(int &index, IntPointData &data, Tensor &Stress){
    
}   