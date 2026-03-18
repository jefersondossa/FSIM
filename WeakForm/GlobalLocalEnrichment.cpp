#include "GlobalLocalEnrichment.h"

GlobalLocalEnrichment::GlobalLocalEnrichment(int matid, int dimension, double young, double poisson) : WeakForm() {
    this->fMatId = matid;
    fNState = dimension;
    fYoungModulus = young;
    fPoissonRatio = poisson;
    fConstitutiveMatrix.resize(3,3);
    fConstitutiveMatrix.setZero();

    double fBulkModulus = fYoungModulus / (2. * (1.-fPoissonRatio));
    double fShearModulus = fYoungModulus / (2. * (1.+fPoissonRatio));
    double k = fYoungModulus / (1. - fPoissonRatio * fPoissonRatio);
    fConstitutiveMatrix(0,0) = k;
    fConstitutiveMatrix(0,1) = k * fPoissonRatio;
    fConstitutiveMatrix(1,0) = k * fPoissonRatio;
    fConstitutiveMatrix(1,1) = k;
    fConstitutiveMatrix(2,2) = k * (1. - fPoissonRatio) * 0.5;  
};

void GlobalLocalEnrichment::ComputeStiffness(int &index, IntPointData &localdata, IntPointData &globaldata, MatrixDouble &Stiffness){
    
    if (!localdata.fNeedsSol){
        localdata.fNeedsSol = true;
        localdata.fSol.resize(fNState);
    }
    if (!globaldata.fNeedsSol){
        globaldata.fNeedsSol = true;
        globaldata.fSol.resize(fNState);
    }

    double WJ = localdata.fWeight * localdata.fJacA0;
    int nphi = globaldata.fPhi.size();

    MatrixDouble matB(3,2*nphi);
    MatrixDouble matBEnr(3,2*nphi);
    MatrixDouble matBTot(3, (matB.cols()+matBEnr.cols()));

    matB.setZero();
    matBEnr.setZero();
    matBTot.setZero();

    double uXInterp = localdata.fSol[0];
    double uYInterp = localdata.fSol[1];
    double dUxdx = localdata.fDSolDx(0,0);
    double dUydy = localdata.fDSolDx(1,1);
    double dUxdy = localdata.fDSolDx(0,1);
    double dUydx = localdata.fDSolDx(1,0);

    for (int j = 0; j < nphi; j++){

        matB(0,2*j  ) = globaldata.fDPhiX0(0,j);
        matB(1,2*j+1) = globaldata.fDPhiX0(1,j);
        matB(2,2*j  ) = globaldata.fDPhiX0(1,j);
        matB(2,2*j+1) = globaldata.fDPhiX0(0,j);

        matBEnr(0,2*j  ) = globaldata.fDPhiX0(0,j)*uXInterp + globaldata.fPhi(j)*dUxdx;
        matBEnr(1,2*j+1) = globaldata.fDPhiX0(1,j)*uYInterp + globaldata.fPhi(j)*dUydy;
        matBEnr(2,2*j  ) = globaldata.fDPhiX0(1,j)*uXInterp + globaldata.fPhi(j)*dUxdy;
        matBEnr(2,2*j+1) = globaldata.fDPhiX0(0,j)*uYInterp + globaldata.fPhi(j)*dUydx;
    }

    //std::cout << "Mat B =\n"<< matB << std::endl;
    //std::cout << "Mat BEnr =\n"<< matBEnr << std::endl;

    matBTot << matB, matBEnr;

    //std::cout << "Mat BTot =\n"<< matBTot.transpose() << std::endl;

    Stiffness += matBTot.transpose() * fConstitutiveMatrix * matBTot * WJ;

    //std::cout << "Stiffness =\n"<< Stiffness << std::endl;
};

void GlobalLocalEnrichment::ComputeResidual(int &index, IntPointData &localdata, IntPointData &globaldata, VecDouble &Rhs){

    double WJ = localdata.fWeight * localdata.fJacA0;
    int nphi = globaldata.fPhi.size();
    
    double nphiL = globaldata.fPhi.size();
    double nphiG = globaldata.fPhi.size();
    VecDouble forcingF(1);
    VecDouble x_ = localdata.fX;
    if (fForceFunction) fForceFunction(x_,forcingF);

    for (int i = 0; i < nphiL; i++){
        Rhs[i] -= WJ * localdata.fPhi[i] * globaldata.fSol[0];
    }
    for (int i = 0; i < nphiG; i++){
        Rhs[nphiL+i] += (forcingF[0] * globaldata.fPhi[i] - globaldata.fPhi[i] * globaldata.fSol[0]) * WJ;
    }
};



int GlobalLocalEnrichment::VariableIndex(const std::string &name) const{
    
    if(!strcmp("GlobalLocalEnrichment",name.c_str()))    return 1;
    
    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int GlobalLocalEnrichment::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
        return fNState;

    default:
        PanicButton();
        return -1;
    }
};

void GlobalLocalEnrichment::Solution(IntPointData &data, int var, VecDouble &Sol) {

    //Solution
    if (var == 1){
        for (int i = 0; i < fNState; i++){
            Sol[i] = data.fSol[i];
        }
        return;
    };

}; 