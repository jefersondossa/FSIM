#include "TransientElasticity2D.h"

TransientElasticity2D::TransientElasticity2D(int matid, double young, double poisson, bool planes, double damp, double dens, double dt) : 
                        Elasticity2D(matid,young,poisson,planes) {
    this->fMatId = matid;
    fDimension = 2;
    fNState = 2;
    fDamping = damp;
    fDensity = dens;
    fTimeStep = dt;
};


void TransientElasticity2D::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){

    //Static Stiffness
    Elasticity2D::ComputeStiffness(index,data,Stiffness);
    
}

void TransientElasticity2D::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    //Static Residual
    Elasticity2D::ComputeResidual(index,data,Rhs);
    
};

void TransientElasticity2D::ComputeError(IntPointData &data, VecDouble &errors){
    PanicButton();
}

int TransientElasticity2D::VariableIndex(const std::string &name) const{
    
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

    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int TransientElasticity2D::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
    case 8:
    case 15:
    case 16:
    case 17:
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
        return 1;

    default:
        PanicButton();
        return -1;
    }
};

void TransientElasticity2D::Solution(IntPointData &data, int var, VecDouble &Sol) {

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
        VecDouble epsilon(3);
        epsilon[0] = data.fDSolDx(0,0);
        epsilon[1] = data.fDSolDx(1,1);
        epsilon[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        auto sigma = fConstitutiveMatrix * epsilon;
        Sol[0] = sigma[0];
        Sol[1] = sigma[1];
        Sol[2] = sigma[2];
        return;
    };

    //Strain
    if (var == 17){
        Sol[0] = data.fDSolDx(0,0);
        Sol[1] = data.fDSolDx(1,1);
        Sol[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);
        return;
    };

}; 

