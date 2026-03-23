#include "Elasticity2D.h"

Elasticity2D::Elasticity2D(int matid, double young, double poisson, bool planes, double thick) : WeakForm() {
    this->fMatId = matid;
    fDimension = 2;
    fNState = 2;
    fYoungModulus = young;
    fPoissonRatio = poisson;
    fPlaneStress = planes;
    fThickness = thick;
    fConstitutiveMatrix.resize(3,3);
    fConstitutiveMatrix.setZero();

    if (fPlaneStress){//Plane Stress Matrix
        double fBulkModulus = fYoungModulus / (2. * (1.-fPoissonRatio));
        double fShearModulus = fYoungModulus / (2. * (1.+fPoissonRatio));
        double k = fYoungModulus / (1. - fPoissonRatio * fPoissonRatio);
        fConstitutiveMatrix(0,0) = k;
        fConstitutiveMatrix(0,1) = k * fPoissonRatio;
        fConstitutiveMatrix(1,0) = k * fPoissonRatio;
        fConstitutiveMatrix(1,1) = k;
        fConstitutiveMatrix(2,2) = k * (1. - fPoissonRatio) * 0.5;
        // double alpha = (3.*fBulkModulus - 2.*fShearModulus) / (3.*fBulkModulus + 4.*fShearModulus);
        // fConstitutiveMatrix(0,0) = fConstitutiveMatrix(1,1) = 1. + alpha;
        // fConstitutiveMatrix(0,1) = fConstitutiveMatrix(1,0) = alpha;
        // fConstitutiveMatrix(2,2) = 0.5;
        // fConstitutiveMatrix *= 2.*fShearModulus;
    } else {//Plane Strain Matrix
        double fBulkModulus = fYoungModulus / (3. * (1.-2.*fPoissonRatio));
        double fShearModulus = fYoungModulus / (2. * (1.+fPoissonRatio));
        double aux = fYoungModulus /(( 1. + fPoissonRatio)*(1.-2.*fPoissonRatio));
        fConstitutiveMatrix(0,0) = (1.-fPoissonRatio) * aux;
        fConstitutiveMatrix(0,1) = aux * fPoissonRatio;
        fConstitutiveMatrix(1,0) = aux * fPoissonRatio;
        fConstitutiveMatrix(1,1) = (1.-fPoissonRatio) * aux;
        fConstitutiveMatrix(2,2) = fShearModulus;
    }
};


void Elasticity2D::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){

    if (!data.fNeedsDSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
        data.fNeedsSol = true;
        data.fSol.resize(fNState);
    }

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index] * fThickness;
    int nphi = data.fPhi.size();
    MatrixDouble matB(3,2*nphi);
    matB.setZero();
    auto matBT=matB.transpose();

    for (int j = 0; j < nphi; j++){
        matB(0,2*j  ) = data.fDPhiX0(0,j);
        matB(1,2*j+1) = data.fDPhiX0(1,j);
        matB(2,2*j  ) = data.fDPhiX0(1,j);
        matB(2,2*j+1) = data.fDPhiX0(0,j);
    }
    
    Stiffness += matB.transpose() * fConstitutiveMatrix * matB * WJ;

    // std::cout << "Stiffness =\n"<< Stiffness << std::endl;
}

void Elasticity2D::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    int nphi = data.fPhi.size();

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index] * fThickness;
    // MatrixDouble matB(3,2*nphi);
    // matB.setZero();
   

    auto force = fForceFunction;
    VecDouble forcingF(fDimension);
    forcingF.setZero();
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);
    
    // for (int j = 0; j < nphi; j++){
    //     matB(0,fDimension*j  ) = data.fDPhiX0(0,j);
    //     matB(1,fDimension*j+1) = data.fDPhiX0(1,j);
    //     matB(2,fDimension*j  ) = data.fDPhiX0(1,j);
    //     matB(2,fDimension*j+1) = data.fDPhiX0(0,j);
    // }
    
    // VecDouble strain(3);
    // strain.setZero();
    // strain[0] = data.fDSolDx(0,0);
    // strain[1] = data.fDSolDx(1,1);
    // strain[2] = (data.fDSolDx(0,1)+data.fDSolDx(1,0));
    // VecDouble stress = fConstitutiveMatrix * strain;

    // Rhs -= matB.transpose() * stress * WJ;

    for (int i = nphi; i--; ){
        double shapeFi = data.fPhi[i];
        //External force
        double Fx = forcingF[0] * shapeFi;
        double Fy = forcingF[1] * shapeFi;
        Rhs[2*i  ] += Fx * WJ;
        Rhs[2*i+1] += Fy * WJ;
    };

    
};

void Elasticity2D::ComputeError(IntPointData &data, VecDouble &errors){
    errors.resize(4);

    VecDouble uExact(fDimension);
    MatrixDouble DuExact(fDimension,fDimension);
    VecDouble x_ = data.fX;
    fExactSol(x_,uExact,DuExact);

    //L2 displacement
    errors[0] += ((uExact[0]-data.fSol[0])*(uExact[0]-data.fSol[0]) + 
                  (uExact[1]-data.fSol[1])*(uExact[1]-data.fSol[1]))
                  * data.fWeight * data.fJacA0 ;

    VecDouble exactStrain(3);
    exactStrain(0) = DuExact(0,0);
    exactStrain(1) = DuExact(1,1);
    exactStrain(2) = 0.5 * (DuExact(1,0) + DuExact(0,1));
    auto exactStress = fConstitutiveMatrix * exactStrain;

    VecDouble StrainMEF(3);
    StrainMEF(0) = data.fDSolDx(0,0);
    StrainMEF(1) = data.fDSolDx(1,1);
    StrainMEF(2) = 0.5 * (data.fDSolDx(1,0) + data.fDSolDx(0,1));
    auto StressMEF = fConstitutiveMatrix * StrainMEF;

    double sigx = StressMEF[0] - exactStress[0];
    double sigy = StressMEF[1] - exactStress[1];
    double sigxy = StressMEF[2] - exactStress[2];

    // Energy norm
    errors[1] = (sigx*(StrainMEF[0]-exactStrain[0])+sigy*(StrainMEF[1]-exactStrain[1])+2.*sigxy*(StrainMEF[2]-exactStrain[2]));
	
	// erro em norma L2 em tensoes
    errors[2] = sigx*sigx + sigy*sigy + 2.*sigxy*sigxy;
    
	// erro estimado na norma H1
    double SemiH1 =0.;
    for(int i = 0; i < 2; i++) for(int j = 0; j < 2; j++) SemiH1 += (data.fDSolDx(i,j) - DuExact(i,j)) * (data.fDSolDx(i,j) - DuExact(i,j));
	errors[3] = errors[0] + SemiH1;
}

int Elasticity2D::VariableIndex(const std::string &name) const{
    
    if(!strcmp("Displacement",name.c_str()))     return 1;
    if(!strcmp("SigmaX",name.c_str()))           return 2;
    if(!strcmp("SigmaY",name.c_str()))           return 3;
    if(!strcmp("TauXY",name.c_str()))            return 4;
    if(!strcmp("EpsilonX",name.c_str()))         return 5;
    if(!strcmp("EpsilonY",name.c_str()))         return 6;
    if(!strcmp("EpsilonXY",name.c_str()))        return 7;
    if(!strcmp("ExactDisplacement",name.c_str()))return 8;
    if(!strcmp("ExactSigmaX",name.c_str()))      return 9;
    if(!strcmp("ExactSigmaY",name.c_str()))      return 10;
    if(!strcmp("ExactTauXY",name.c_str()))       return 11;
    if(!strcmp("ExactEpsilonX",name.c_str()))    return 12;
    if(!strcmp("ExactEpsilonY",name.c_str()))    return 13;
    if(!strcmp("ExactEpsilonXY",name.c_str()))   return 14;
    if(!strcmp("ExactForce",name.c_str()))       return 15;
    if(!strcmp("Stress",name.c_str()))           return 16;
    if(!strcmp("Strain",name.c_str()))           return 17;
    if(!strcmp("SigmaZ",name.c_str()))           return 18;
    if(!strcmp("DeltaStrain",name.c_str()))      return 19;
    if(!strcmp("Pressure",name.c_str()))         return 20;
    if(!strcmp("J2",name.c_str()))               return 21;
    if(!strcmp("Compliance",name.c_str()))       return 22;

    // std::cout << "Post Process variable not implemented \n";
    // PanicButton();
    return -1;
};

int Elasticity2D::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
    case 8:
    case 15:
    case 16:
    case 17:
    case 19:
        return 3;
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
    case 18:
    case 20:
    case 21:
    case 22:
        return 1;

    default:
        // PanicButton();
        return -1;
    }
};

void Elasticity2D::Solution(IntPointData &data, int var, VecDouble &Sol) {

    //Displacement
    if (var == 1){
        Sol[0] = data.fSol[0];
        Sol[1] = data.fSol[1];
        Sol[2] = 0.;
        return;
    };

    //Sigma X
    if (var == 2){
        VecDouble epsilon(3);
        epsilon[0] = data.fDSolDx(0,0);
        epsilon[1] = data.fDSolDx(1,1);
        epsilon[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        if (fPlaneStress){
            double k = fYoungModulus / (1.-fPoissonRatio*fPoissonRatio);
            Sol[0] = k * (epsilon[0] + fPoissonRatio * epsilon[1]);
        } else {
            double k = fYoungModulus / ((1.+fPoissonRatio)*(1.-2.*fPoissonRatio));
            Sol[0] = k * ((1.-fPoissonRatio) * epsilon[0] + fPoissonRatio * epsilon[1]);
        }        
        return;
    };

    //Sigma Y
    if (var == 3){
        VecDouble epsilon(3);
        epsilon[0] = data.fDSolDx(0,0);
        epsilon[1] = data.fDSolDx(1,1);
        epsilon[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        if (fPlaneStress){
            double k = fYoungModulus / (1.-fPoissonRatio*fPoissonRatio);
            Sol[0] = k * (fPoissonRatio * epsilon[0] + epsilon[1]);
        } else {
            double k = fYoungModulus / ((1.+fPoissonRatio)*(1.-2.*fPoissonRatio));
            Sol[0] = k * (fPoissonRatio * epsilon[0] + (1.-fPoissonRatio) * epsilon[1]);
        }      
        return;
    };

    //Tau XY
    if (var == 4){
        VecDouble epsilon(3);
        epsilon[0] = data.fDSolDx(0,0);
        epsilon[1] = data.fDSolDx(1,1);
        epsilon[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        if (fPlaneStress){
            double k = fYoungModulus / (1.-fPoissonRatio*fPoissonRatio);
            Sol[0] = k * (1.-fPoissonRatio) * epsilon[2];
        } else {
            double k = fYoungModulus / ((1.+fPoissonRatio)*(1.-2.*fPoissonRatio));
            Sol[0] = k * (1.-2.*fPoissonRatio) * epsilon[2];
        }      
        return;
    };

    //Epsilon X
    if (var == 5){
        Sol[0] = data.fDSolDx(0,0);
        return;
    };
    //Epsilon Y
    if (var == 6){
        Sol[0] = data.fDSolDx(1,1);
        return;
    };
    //Epsilon XY
    if (var == 7){
        Sol[0] = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        return;
    };


    VecDouble forcingF(fDimension);
    VecDouble x_ = data.fX;
    if (fForceFunction) fForceFunction(x_,forcingF);

    VecDouble disp(fDimension);
    MatrixDouble gradDisp(fDimension,3);
    if (fExactSol) fExactSol(x_,disp,gradDisp);

    //Exact Displacement
    if (var == 8){
        Sol[0] = disp[0];
        Sol[1] = disp[1];
        Sol[2] = 0.;
        return;
    };
    
    //Exact Sigma X
    if (var == 9){
        VecDouble epsilon(3);
        epsilon[0] = gradDisp(0,0);
        epsilon[1] = gradDisp(1,1);
        epsilon[2] = gradDisp(0,1)+gradDisp(1,0);
        if (fPlaneStress){
            double k = fYoungModulus / (1.-fPoissonRatio*fPoissonRatio);
            Sol[0] = k * (epsilon[0] + fPoissonRatio * epsilon[1]);
        } else {
            double k = fYoungModulus / ((1.+fPoissonRatio)*(1.-2.*fPoissonRatio));
            Sol[0] = k * ((1.-fPoissonRatio) * epsilon[0] + fPoissonRatio * epsilon[1]);
        }        
        return;
    };

    //Exact Sigma Y
    if (var == 10){
        VecDouble epsilon(3);
        epsilon[0] = gradDisp(0,0);
        epsilon[1] = gradDisp(1,1);
        epsilon[2] = gradDisp(0,1)+gradDisp(1,0);
        if (fPlaneStress){
            double k = fYoungModulus / (1.-fPoissonRatio*fPoissonRatio);
            Sol[0] = k * (fPoissonRatio * epsilon[0] + epsilon[1]);
        } else {
            double k = fYoungModulus / ((1.+fPoissonRatio)*(1.-2.*fPoissonRatio));
            Sol[0] = k * (fPoissonRatio * epsilon[0] + (1.-fPoissonRatio) * epsilon[1]);
        }      
        return;
    };

    //Exact Tau XY
    if (var == 11){
        VecDouble epsilon(3);
        epsilon[0] = gradDisp(0,0);
        epsilon[1] = gradDisp(1,1);
        epsilon[2] = gradDisp(0,1)+gradDisp(1,0);
        if (fPlaneStress){
            double k = fYoungModulus / (1.-fPoissonRatio*fPoissonRatio);
            Sol[2] = k * (1.-fPoissonRatio) * epsilon[2];
        } else {
            double k = fYoungModulus / ((1.+fPoissonRatio)*(1.-2.*fPoissonRatio));
            Sol[2] = k * (1.-2.*fPoissonRatio) * epsilon[2];
        }      
        return;
    };

    //Exact Epsilon X
    if (var == 12){
        Sol[0] = gradDisp(0,0);
        return;
    };
    //Exact Epsilon Y
    if (var == 13){
        Sol[0] = gradDisp(1,1);
        return;
    };
    //Exact Epsilon XY
    if (var == 14){
        Sol[0] = gradDisp(0,1)+gradDisp(1,0);
        return;
    };

    //Exact Force
    if (var == 15){
        Sol[0] = forcingF[0];
        Sol[1] = forcingF[1];
        Sol[2] = 0.;
        return;
    };

    //Stress
    if (var == 16){
        VecDouble epsilon(3);
        epsilon[0] = data.fDSolDx(0,0);
        epsilon[1] = data.fDSolDx(1,1);
        epsilon[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        if (fPlaneStress){
            double fBulkModulus = fYoungModulus / (2. * (1.-fPoissonRatio));
            double fShearModulus = fYoungModulus / (2. * (1.+fPoissonRatio));
            double alpha = (3.*fBulkModulus - 2.*fShearModulus) / (3.*fBulkModulus + 4.*fShearModulus);
            MatrixDouble MatAux(3,3);
            MatAux.setZero();
            MatAux(0,0) = MatAux(1,1) = 1. + alpha;
            MatAux(0,1) = MatAux(1,0) = alpha;
            MatAux(2,2) = 0.5;
            MatAux *= 2.*fShearModulus;
            Sol = MatAux*epsilon;
            // double k = fYoungModulus / (1.-fPoissonRatio*fPoissonRatio);
            // Sol[0] = k * (epsilon[0] + fPoissonRatio * epsilon[1]);
            // Sol[1] = k * (fPoissonRatio * epsilon[0] + epsilon[1]);
            // Sol[2] = k * (1.-fPoissonRatio) * epsilon[2];
        } else {
            double G = fYoungModulus / (2. * ( 1. + fPoissonRatio));
            
            double k = fYoungModulus / ((1.+fPoissonRatio)*(1.-2.*fPoissonRatio));
            Sol[0] = k * ((1.-fPoissonRatio) * epsilon[0] + fPoissonRatio * epsilon[1]);
            Sol[1] = k * (fPoissonRatio * epsilon[0] + (1.-fPoissonRatio) * epsilon[1]);
            Sol[2] = G * epsilon[2];
        }
        return;
    };

    //Strain
    if (var == 17){
        Sol[0] = data.fDSolDx(0,0);
        Sol[1] = data.fDSolDx(1,1);
        Sol[2] = (data.fDSolDx(0,1)+data.fDSolDx(1,0));
        return;
    };

    //Stress Z
    if (var == 18){
        VecDouble epsilon(3);
        epsilon[0] = data.fDSolDx(0,0);
        epsilon[1] = data.fDSolDx(1,1);
        epsilon[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        if (fPlaneStress){
            Sol[0] = 0.;
        } else {
            double k = fYoungModulus / ((1.+fPoissonRatio)*(1.-2.*fPoissonRatio));
            Sol[0] = k * (fPoissonRatio * epsilon[0] + fPoissonRatio * epsilon[1]);
        }      
        return;
    };

    //Delta Strain
    if (var == 19){
        Sol[0] = data.fDSolDx(0,0) - data.fDSolDxPrev(0,0);
        Sol[1] = data.fDSolDx(1,1) - data.fDSolDxPrev(1,1);
        Sol[2] = ((data.fDSolDx(0,1)+data.fDSolDx(1,0))-(data.fDSolDxPrev(0,1)+data.fDSolDxPrev(1,0)));
        // int index = data.fIndex;
        // Sol[0] = data.fDSolDx(0,0) - data.fElasticStrain[index].fXX();
        // Sol[1] = data.fDSolDx(1,1) - data.fElasticStrain[index].fYY();
        // Sol[2] = ((data.fDSolDx(0,1)+data.fDSolDx(1,0))+(data.fElasticStrain[index].fXY()*2.));
        // Sol[0] = data.fDSolDxPrev(0,0) + data.fDSolDx(0,0);
        // Sol[1] = data.fDSolDxPrev(1,1) + data.fDSolDx(1,1);
        // Sol[2] = (data.fDSolDxPrev(0,1)+data.fDSolDxPrev(1,0))+(data.fDSolDx(0,1)+data.fDSolDx(1,0));
        return;
    };

    //Pressure
    if (var == 20){
        VecDouble epsilon(3);
        epsilon[0] = data.fDSolDx(0,0);
        epsilon[1] = data.fDSolDx(1,1);
        epsilon[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        if (fPlaneStress){
            PanicButton();
        } else {
            double G = fYoungModulus / (2. * ( 1. + fPoissonRatio));
            double k = fYoungModulus / ((1.+fPoissonRatio)*(1.-2.*fPoissonRatio));
            Tensor3D Stress;
            Stress.fXX() = k * ((1.-fPoissonRatio) * epsilon[0] + fPoissonRatio * epsilon[1]);
            Stress.fYY() = k * (fPoissonRatio * epsilon[0] + (1.-fPoissonRatio) * epsilon[1]);
            Stress.fZZ() = k * (fPoissonRatio * epsilon[0] + fPoissonRatio * epsilon[1]);
            Stress.fXY() = G * epsilon[2];
            Sol[0] = Stress.Trace()/3.;
        }
        return;
    };
    
    //J2
    if (var == 21){
        VecDouble epsilon(3);
        epsilon[0] = data.fDSolDx(0,0);
        epsilon[1] = data.fDSolDx(1,1);
        epsilon[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        if (fPlaneStress){
            PanicButton();
        } else {
            double G = fYoungModulus / (2. * ( 1. + fPoissonRatio));
            double k = fYoungModulus / ((1.+fPoissonRatio)*(1.-2.*fPoissonRatio));
            Tensor3D Stress;
            Stress.fXX() = k * ((1.-fPoissonRatio) * epsilon[0] + fPoissonRatio * epsilon[1]);
            Stress.fYY() = k * (fPoissonRatio * epsilon[0] + (1.-fPoissonRatio) * epsilon[1]);
            Stress.fZZ() = k * (fPoissonRatio * epsilon[0] + fPoissonRatio * epsilon[1]);
            Stress.fXY() = G * epsilon[2];
            Sol[0] = Stress.J2();
        }
        return;
    };

    //Compliance
    if (var == 22){
        if(!data.fStiffnessMatrix.has_value()){
            Sol[0] = std::numeric_limits<float>::max();
            PanicButton();
            return;
        }
        const auto compliance = data.fSolNodes.transpose() * (*data.fStiffnessMatrix) * data.fSolNodes;
        Sol[0] = compliance(0, 0);
        return;
    }
};


MatrixDouble &Elasticity2D::ConstitutiveMatrix(){
    return fConstitutiveMatrix;
}