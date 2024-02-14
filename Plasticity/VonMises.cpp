#include "VonMises.h"

VonMises::VonMises(WeakForm *elast) : PlasticityModel(elast){
    
}

void VonMises::ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor &Stress){
    
    VecDouble deviatory(fNStressComponents);
    MatrixDouble fTangentTensor(fNStressComponents,fNStressComponents);
    fTangentTensor.setZero();
    deviatory.setZero();
    auto devAux = Stress.Deviatory();

    if (fRealDimension == 2){   
        deviatory[0] = devAux.fData[0]; 
        deviatory[1] = devAux.fData[1]; 
        deviatory[2] = devAux.fData[3]; 
    } else if (fRealDimension == 3){
        PanicButton();
    } else {
        PanicButton();
    }
    double vonMisesStress = Stress.DeviatoryNorm();
    double Qtrial = sqrt(1.5) * vonMisesStress;

    double Afactor = 2. * fShearModulus *(1.-3.*fShearModulus*data.fPlasticMultiplier/Qtrial);
    double Bfactor = 6. * fShearModulus * fShearModulus * (data.fPlasticMultiplier/Qtrial - 1./(3.*fShearModulus+fHardening))/(vonMisesStress*vonMisesStress);

    for (int i = 0; i < fNStressComponents; i++){
        for (int j = 0; j < fNStressComponents; j++){
            fTangentTensor(i,j) = Afactor*fIdentity4Dev(i,j) + Bfactor*deviatory[i]*deviatory[j] + fBulkModulus*fIdentity2[i]*fIdentity2[j];  
        }
    }

    // fElasticModel->ConstitutiveMatrix() = fTangentTensor;
    fElasticModel->ComputeStiffness(index,data,Stiffness);

};
    
void VonMises::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs, Tensor &Stress){
    
    // if (data.fPlasticStrain.norm() > 0){
        
        int nphi = data.fPhi.size();

        double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
        MatrixDouble matB(3,2*nphi);
        matB.setZero();

        auto force = fForceFunction;
        VecDouble forcingF(fDimension);
        forcingF.setZero();
        VecDouble x_ = data.fX;
        if (force) force(x_,forcingF);
        
        // for (int j = 0; j < nphi; j++){
        //     matB(0,fDimension*j  ) = data.fDPhiX0(j,0);
        //     matB(1,fDimension*j+1) = data.fDPhiX0(j,1);
        //     matB(2,fDimension*j  ) = data.fDPhiX0(j,1);
        //     matB(2,fDimension*j+1) = data.fDPhiX0(j,0);
        // }
        
        // MatrixDouble elasticStrain(2,2);
        // elasticStrain(0,0) = data.fDSolDx(0,0)+data.fPlasticStrain[index];
        // elasticStrain(1,1) = data.fDSolDx(1,1)+data.fPlasticStrain[index];
        // elasticStrain(0,1) = elasticStrain(1,0) = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        // Tensor EStrain(elasticStrain);
        // double volStrain = EStrain.Trace();
        // double pressure = fBulkModulus * volStrain;
        // auto elasticDev = EStrain.Deviatory();

        // VecDouble strain(3),stress(3);
        // strain[0] = elasticDev.fData[0] / (2. * fShearModulus) + volStrain/2.;
        // strain[1] = elasticDev.fData[1] / (2. * fShearModulus) + volStrain/2.;
        // strain[2] = elasticDev.fData[3] / (fShearModulus);

        // double gamma = data.fPlasticMultiplier[index];
        // data.fPlasticStrain[index] += gamma;

        // double vonMisesStress = Stress.DeviatoryNorm();
        // double Q = sqrt(1.5) * vonMisesStress;
        // double Qtrial = Q + 3.*fShearModulus*gamma;

        // double factor = 2. * fShearModulus * (1. - 3.*fShearModulus*gamma/Qtrial);
        // stress[0] = factor*elasticDev.fData[0] + pressure;
        // stress[1] = factor*elasticDev.fData[1] + pressure;
        // stress[2] = factor*elasticDev.fData[3];
        // stress.setZero();



        // double k = fYoungModulus / (1. - fPoissonRatio * fPoissonRatio);
        // fElasticModel->ConstitutiveMatrix().setZero();
        // fElasticModel->ConstitutiveMatrix()(0,0) = k;
        // fElasticModel->ConstitutiveMatrix()(0,1) = k * fPoissonRatio;
        // fElasticModel->ConstitutiveMatrix()(1,0) = k * fPoissonRatio;
        // fElasticModel->ConstitutiveMatrix()(1,1) = k;
        // fElasticModel->ConstitutiveMatrix()(2,2) = k * (1. - fPoissonRatio) * 0.5;

        // VecDouble strain(3);
        // strain.setZero();
        // strain[0] = data.fDSolDx(0,0);
        // strain[1] = data.fDSolDx(1,1);
        // strain[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0); 

        VecDouble stress(3);// = fElasticModel->ConstitutiveMatrix() * strain;
        // double vonMisesStress = sqrt(2.0 / 3.0) * stress.norm();
        // double deltaEp = (vonMisesStress - 0.45) / (2.0 * fYoungModulus + fHardening);
        // stress -= 2.0 * fYoungModulus * deltaEp * strain / vonMisesStress;
        // data.fPlasticStrain[index] += deltaEp;

        // TensorToVoigt(Stress,stress);
        stress[0] = Stress.fData[0];
        stress[1] = Stress.fData[1];
        stress[2] = Stress.fData[3];

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
       
        // fElasticModel->ComputeResidual(index,data,Rhs);
    // }else{
        // fElasticModel->ComputeResidual(index,data,Rhs);
    // }
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
        Sol[0] = 1.;
        return;
    };
};
double VonMises::YieldFunction(int &index, IntPointData &data, Tensor &Stress){
    //Eq (7.74)
    double plasticstrain = data.fPlasticStrain[index];
    double fYield;
    fUniaxialYield(plasticstrain,fYield,fHardening);
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

double VonMises::PlasticMultiplier(int &index, IntPointData &data, Tensor &Stress){
    //Newton-Raphson to find plastic multiplier

    //Plastic Multiplier - Box 7.4 Souza Neto
    double dGamma = 0.;
    double Qtrial = sqrt(1.5) * Stress.DeviatoryNorm();
    double sigmay = 0.;
    fUniaxialYield(data.fPlasticStrain[index],sigmay,fHardening);
    double deltaGamma = 0.;
    double PhiTil = Qtrial-sigmay;
    
    while (fabs(PhiTil/sigmay) > 1.e-5){
        double d = -3. * fShearModulus - fHardening;
        dGamma -= PhiTil/d;
        fUniaxialYield(data.fPlasticStrain[index]+dGamma,sigmay,fHardening);
        PhiTil = Qtrial - 3.*fShearModulus*dGamma-sigmay;
    }
    
    return dGamma;
}

void VonMises::UpdateStateVariables(int &index, IntPointData &data, Tensor &Stress){

    auto dev = Stress.Deviatory();
    double Qtrial = sqrt(1.5) * Stress.DeviatoryNorm();
    auto hydrostatic = Stress.Hydrostatic();
    dev *= (1. - data.fPlasticMultiplier * 3 * fShearModulus / Qtrial);
    Stress = dev + hydrostatic;
    data.fPlasticStrain[index] += data.fPlasticMultiplier;
}   