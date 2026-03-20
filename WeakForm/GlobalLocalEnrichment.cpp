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

    // int nphi = globaldata.fPhi.size();
    // double WJ = localdata.fWeight * localdata.fJacA0;

    // MatrixDouble matB(3,2*nphi);
    // MatrixDouble matBEnr(3,2*nphi);
    // MatrixDouble matBTot(3, (matB.cols()+matBEnr.cols()));
    
    // matB.setZero();
    // matBEnr.setZero();
    // matBTot.setZero();

    // auto force = fForceFunction;
    // VecDouble forcingF(fDimension);
    // forcingF.setZero();
    // VecDouble x_ = globaldata.fX;
    // if (force) force(x_,forcingF);

    // double uXInterp = localdata.fSol[0];
    // double uYInterp = localdata.fSol[1];
    // double dUxdx = localdata.fDSolDx(0,0);
    // double dUydy = localdata.fDSolDx(1,1);
    // double dUxdy = localdata.fDSolDx(0,1);
    // double dUydx = localdata.fDSolDx(1,0);
    
    // for (int j = 0; j < nphi; j++){
    //     matB(0,2*j  ) = globaldata.fDPhiX0(0,j);
    //     matB(1,2*j+1) = globaldata.fDPhiX0(1,j);
    //     matB(2,2*j  ) = globaldata.fDPhiX0(1,j);
    //     matB(2,2*j+1) = globaldata.fDPhiX0(0,j);

    //     matBEnr(0,2*j  ) = globaldata.fDPhiX0(0,j)*uXInterp + globaldata.fPhi(j)*dUxdx;
    //     matBEnr(1,2*j+1) = globaldata.fDPhiX0(1,j)*uYInterp + globaldata.fPhi(j)*dUydy;
    //     matBEnr(2,2*j  ) = globaldata.fDPhiX0(1,j)*uXInterp + globaldata.fPhi(j)*dUxdy;
    //     matBEnr(2,2*j+1) = globaldata.fDPhiX0(0,j)*uYInterp + globaldata.fPhi(j)*dUydx;
    // }

    // matBTot << matB, matBEnr;
    
    // VecDouble strain(3);
    // strain.setZero();
    // strain[0] = globaldata.fDSolDx(0,0);
    // strain[1] = globaldata.fDSolDx(1,1);
    // strain[2] = (globaldata.fDSolDx(0,1)+globaldata.fDSolDx(1,0));
    // VecDouble stress = fConstitutiveMatrix * strain;

    // Rhs -= matBTot.transpose() * stress * WJ;

    // for (int i = nphi; i--; ){
    //     double shapeFi = globaldata.fPhi[i];
    //     //External force
    //     double Fx = uXInterp * shapeFi; 
    //     double Fy = uYInterp * shapeFi;
    //     Rhs[2*nphi + 2*i] += Fx * WJ;
    //     Rhs[2*nphi + 2*i+1] += Fy * WJ;
    // };
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