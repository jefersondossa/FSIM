#include "ElasticityPositional2D.h"

ElasticityPositional2D::ElasticityPositional2D(int matid, REAL young, REAL poisson, REAL thick, bool planes) : WeakForm() {
    this->fMatId = matid;
    fDimension = 2;
    fNState = 2;
    fYoungModulus = young;
    fPoissonRatio = poisson;
    fPlaneStress = planes;
    this->fType = WeakFormType::kElasticityPositional2D;
    fConstitutiveMatrix.resize(3,3);
    fConstitutiveMatrix.setZero();
    fThickness = thick;

    if (fPlaneStress){
        fConstitutiveMatrix(0,0) = fYoungModulus / (1.0-(fPoissonRatio*fPoissonRatio));
        fConstitutiveMatrix(0,1) = fYoungModulus / (1.0-(fPoissonRatio*fPoissonRatio))*fPoissonRatio;
        fConstitutiveMatrix(1,1) = fYoungModulus / (1.0-(fPoissonRatio*fPoissonRatio));
        fConstitutiveMatrix(1,0) = fYoungModulus / (1.0-(fPoissonRatio*fPoissonRatio))*fPoissonRatio;
        fConstitutiveMatrix(2,2) = 2.0 * (fYoungModulus / (2.0 * (1.0+fPoissonRatio)));
    } else {
        const REAL prop1 = fYoungModulus / ((1.0 + fPoissonRatio) * (1.0 - 2.0 * fPoissonRatio));
        const REAL prop2 = 1.0 - fPoissonRatio;
        const REAL prop3 = fYoungModulus / (1.0 + fPoissonRatio);
        fConstitutiveMatrix(0,0) = prop1*prop2;
        fConstitutiveMatrix(0,1) = prop1*fPoissonRatio;
        fConstitutiveMatrix(1,1) = prop1*prop2;
        fConstitutiveMatrix(1,0) = prop1*fPoissonRatio;
        fConstitutiveMatrix(2,2) = prop3;
    }
    SPKStress.resize(2,2);
    SPKStress.setZero();
    
};

void ElasticityPositional2D::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    if (!data.fNeedsDSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
    }
    int nphi = data.fPhi.size();
    REAL WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index] * fThickness;
    
    auto dphi_dx = data.fDPhiX0.transpose();

    //COMPUTE A1
    MatrixDouble dx_dxsi = data.fA0;
    MatrixDouble dy_dxsi = data.fA1;
    REAL j0 = data.fJacA0;

    //dy_dx
    MatrixDouble dy_dx = dy_dxsi * data.fA0Inv;
    //jacobian
    REAL jac = dy_dx.determinant();

    //Green-Lagrange strain tensor
    MatrixDouble E(fDimension,fDimension);
    for (int i = 0; i < fDimension; i++)
        for (int j = 0; j < fDimension; j++)
            E(i,j) = 0.5 * (dy_dx(0,i) * dy_dx(0,j) + dy_dx(1,i) * dy_dx(1,j));
    E(0,0) -= 0.5; E(1,1) -= 0.5;

    VecDouble EAux(3);
    EAux[0] = E(0,0);
    EAux[1] = E(1,1);
    EAux[2] = E(0,1);

    VecDouble stress = fConstitutiveMatrix * EAux;
    SPKStress(0,0) = stress[0];
    SPKStress(1,1) = stress[1];
    SPKStress(0,1) = SPKStress(1,0) = stress[2];

    //element rhs vector
    for (int a = 0; a < nphi; a++){
        for (int k = 0; k < fDimension; k++){
            MatrixDouble dE_dyak(fDimension,fDimension);
            dE_dyak.setZero();
            for (int i = 0; i < fDimension; i++)
                for (int j = 0; j < fDimension; j++)
                    dE_dyak(i,j) = 0.5 * (dphi_dx(a,i) * dy_dx(k,j) + dy_dx(k,i) * dphi_dx(a,j));

            //element tangent matrix
            for (int b = 0; b < nphi; b++){
                for (int l = 0; l < fDimension; l++){
                    MatrixDouble dE_dybl(fDimension,fDimension);
                    dE_dybl.setZero();
                    for (int i = 0; i < fDimension; i++)
                        for (int j = 0; j < fDimension; j++)
                            dE_dybl(i,j) = 0.5 * (dphi_dx(b,i) * dy_dx(l,j) + dy_dx(l,i) * dphi_dx(b,j));

                    MatrixDouble d2E_dyakbl(fDimension,fDimension);
                    d2E_dyakbl.setZero();
                    if (k==l)
                        for (int i = 0; i < fDimension; i++)
                            for (int j = 0; j < fDimension; j++)
                                d2E_dyakbl(i,j) = 0.5 * (dphi_dx(a,i) * dphi_dx(b,j) + dphi_dx(b,i) * dphi_dx(a,j));
                    else
                        for (int i = 0; i < fDimension; i++)
                            for (int j = 0; j < fDimension; j++)
                                d2E_dyakbl(i,j) = 0.0; 

                    MatrixDouble dS_dybl(fDimension,fDimension);
                    VecDouble auxdedyb(3);
                    auxdedyb[0] = dE_dybl(0,0);
                    auxdedyb[1] = dE_dybl(1,1);
                    auxdedyb[2] = dE_dybl(0,1);
                    stress = fConstitutiveMatrix * auxdedyb;
                    dS_dybl(0,0) = stress[0];
                    dS_dybl(1,1) = stress[1];
                    dS_dybl(0,1) = dS_dybl(1,0) = stress[2];

                    //elastic and geometric componentes of tangent matrix
                    REAL e = 0.0;
                    for (int i = 0; i < fDimension; i++)
                        for (int j = 0; j < fDimension; j++)
                            e += dS_dybl(i,j) * dE_dyak(i,j) + SPKStress(i,j) * d2E_dyakbl(i,j);

                    //mass matrix
                    REAL m = 0.;
                    // (k==l)? m = (1.0 / (0.25 * dTime_)) *0* data.fPhi[a] * data.fPhi[b] : m = 0.0;

                    Stiffness(2 * a + k,2 * b + l) += (e+m) * j0 * data.fWeight * fThickness;
                }
            }
        }
    }


}

void ElasticityPositional2D::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    auto dphi_dx = data.fDPhiX0.transpose();
    int nphi = data.fPhi.size();
    //COMPUTE A1
    auto dx_dxsi = data.fA0;
    auto dy_dxsi = data.fA1;
    auto j0 = data.fJacA0;

    //dy_dx
    auto dy_dx = dy_dxsi * data.fA0Inv;
    //jacobian
    auto jac = data.fJacA1;

    //element rhs vector
    for (int a = 0; a < nphi; a++){
        for (int k = 0; k < fDimension; k++){
            MatrixDouble dE_dyak(fDimension,fDimension);
            dE_dyak.setZero();
            for (int i = 0; i < fDimension; i++)
                for (int j = 0; j < fDimension; j++)
                    dE_dyak(i,j) = 0.5 * (dphi_dx(a,i) * dy_dx(k,j) + dy_dx(k,i) * dphi_dx(a,j));
            
            //internal force
            REAL f = 0.0;
            for (int i = 0; i < fDimension; i++)
                for (int j = 0; j < fDimension; j++)
                    f += SPKStress(i,j) * dE_dyak(i,j);

            Rhs[2 * a + k] -= f * data.fWeight * j0 * fThickness;
        }
    }
};


void ElasticityPositional2D::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs, Tensor3D &Stress){

    auto dphi_dx = data.fDPhiX0;
    int nphi = data.fPhi.size();
    //COMPUTE A1
    auto dx_dxsi = data.fA0;
    auto dy_dxsi = data.fA1;
    auto j0 = data.fJacA0;

    //dy_dx
    auto dy_dx = dy_dxsi * data.fA0Inv;
    //jacobian
    auto jac = data.fJacA1;

    MatrixDouble matStress = Stress.MatrixForm();

    //element rhs vector
    for (int a = 0; a < nphi; a++){
        for (int k = 0; k < fDimension; k++){
            MatrixDouble dE_dyak(fDimension,fDimension);
            dE_dyak.setZero();
            for (int i = 0; i < fDimension; i++)
                for (int j = 0; j < fDimension; j++)
                    dE_dyak(i,j) = 0.5 * (dphi_dx(a,i) * dy_dx(k,j) + dy_dx(k,i) * dphi_dx(a,j));
            
            //internal force
            REAL f = 0.0;
            for (int i = 0; i < fDimension; i++)
                for (int j = 0; j < fDimension; j++)
                    f += matStress(i,j) * dE_dyak(i,j);

            Rhs[2 * a + k] -= f * data.fWeight * j0;
        }
    }
};

void ElasticityPositional2D::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();

}

int ElasticityPositional2D::VariableIndex(const std::string &name) const{
    
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
    if(!strcmp("Material",name.c_str()))         return 100;
    if(!strcmp("DeltaStrain",name.c_str()))      return 19;

    // std::cout << "Post Process variable not implemented \n";
    // PanicButton();
    return -1;
};

int ElasticityPositional2D::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
    case 8:
    case 15:
    case 16:
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
    case 100:
        return 1;

    default:
        // PanicButton();
        return -1;
    }
};

void ElasticityPositional2D::Solution(IntPointData &data, int var, VecDouble &Sol) {

    if (var == 100){
        Sol[0] = fMatId;
        return;
    };

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
        auto sigma = fConstitutiveMatrix * epsilon;
        Sol[0] = sigma[0];
        return;
    };

    //Sigma Y
    if (var == 3){
        VecDouble epsilon(3);
        epsilon[0] = data.fDSolDx(0,0);
        epsilon[1] = data.fDSolDx(1,1);
        epsilon[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        auto sigma = fConstitutiveMatrix * epsilon;
        Sol[0] = sigma[1];
        return;
    };

    //Tau XY
    if (var == 4){
        VecDouble epsilon(3);
        epsilon[0] = data.fDSolDx(0,0);
        epsilon[1] = data.fDSolDx(1,1);
        epsilon[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        auto sigma = fConstitutiveMatrix * epsilon;
        Sol[0] = sigma[2];
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
        auto sigma = fConstitutiveMatrix * epsilon;
        Sol[0] = sigma[0];
        return;
    };

    //Exact Sigma Y
    if (var == 10){
        VecDouble epsilon(3);
        epsilon[0] = gradDisp(0,0);
        epsilon[1] = gradDisp(1,1);
        epsilon[2] = gradDisp(0,1)+gradDisp(1,0);
        auto sigma = fConstitutiveMatrix * epsilon;
        Sol[0] = sigma[1];
        return;
    };

    //Exact Tau XY
    if (var == 11){
        VecDouble epsilon(3);
        epsilon[0] = gradDisp(0,0);
        epsilon[1] = gradDisp(1,1);
        epsilon[2] = gradDisp(0,1)+gradDisp(1,0);
        auto sigma = fConstitutiveMatrix * epsilon;
        Sol[0] = sigma[2];
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
        //dy_dx
        MatrixDouble dy_dx = data.fA1 * data.fA0Inv;
        //jacobian
        REAL jac = dy_dx.determinant();

        //Green-Lagrange strain tensor
        MatrixDouble E(fDimension,fDimension);
        for (int i = 0; i < fDimension; i++)
            for (int j = 0; j < fDimension; j++)
                E(i,j) = 0.5 * (dy_dx(0,i) * dy_dx(0,j) + dy_dx(1,i) * dy_dx(1,j));
        E(0,0) -= 0.5; E(1,1) -= 0.5;
        MatrixDouble constitutive(2,2);
        //Second Piola-Kirchhoff stress tensor
        if (fPlaneStress){
            constitutive(0,0) = fYoungModulus / (1.0-(fPoissonRatio*fPoissonRatio)) * (E(0,0) + fPoissonRatio * E(1,1));
            constitutive(0,1) = 2.0 * (fYoungModulus / (2.0 * (1.0+fPoissonRatio))) * E(0,1);
            constitutive(1,0) = constitutive(0,1);
            constitutive(1,1) = fYoungModulus / (1.0-(fPoissonRatio*fPoissonRatio)) * (E(1,1) + fPoissonRatio * E(0,0));
        } else {
            const REAL prop1 = fYoungModulus / ((1.0 + fPoissonRatio) * (1.0 - 2.0 * fPoissonRatio));
            const REAL prop2 = 1.0 - fPoissonRatio;
            const REAL prop3 = fYoungModulus / (1.0 + fPoissonRatio);
            constitutive(0,0) = prop1 *(prop2 * E(0,0) + fPoissonRatio * E(1,1));    
            constitutive(1,1) = prop1 *(prop2 * E(1,1) + fPoissonRatio * E(0,0));    
            constitutive(0,1) = constitutive(1,0) = prop3 * E(1,0);      
        }
        auto sigma = dy_dx * constitutive * dy_dx.transpose() / jac;
        Sol[0] = sigma(0,0);
        Sol[1] = sigma(1,1);
        Sol[2] = sigma(0,1);
        return;
    };

    if (var == 19){
        //dy_dx
        MatrixDouble dy_dx = data.fA1 * data.fA0Inv;
        MatrixDouble dy_dxPrev = data.fA1Prev * data.fA0Inv;

        //Green-Lagrange strain tensor
        MatrixDouble E(fDimension,fDimension);
        MatrixDouble EPrev(fDimension,fDimension);
        for (int i = 0; i < fDimension; i++){
            for (int j = 0; j < fDimension; j++){
                E(i,j) = 0.5 * (dy_dx(0,i) * dy_dx(0,j) + dy_dx(1,i) * dy_dx(1,j));
                EPrev(i,j) = 0.5 * (dy_dxPrev(0,i) * dy_dxPrev(0,j) + dy_dxPrev(1,i) * dy_dxPrev(1,j));
            }
        }
        E(0,0) -= 0.5; E(1,1) -= 0.5;
        EPrev(0,0) -= 0.5; EPrev(1,1) -= 0.5;

        Sol[0] = E(0,0)-EPrev(0,0);
        Sol[1] = E(1,1)-EPrev(1,1);
        Sol[2] = E(0,1)-EPrev(0,1);
        
        return;
    };


}; 

MatrixDouble &ElasticityPositional2D::ConstitutiveMatrix(){
    return fConstitutiveMatrix;
}