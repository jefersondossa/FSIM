#include "MixedGlobalLocalEnrichment.h"

MixedGlobalLocalEnrichment::MixedGlobalLocalEnrichment(int matid, int dimension, double young, double poisson) : WeakForm() {
    this->fMatId = matid;
    this->fNState = 2;
    fDimension = dimension;
    fYoungModulus = young;
    fPoissonRatio = poisson;
    fBulkModulus = fYoungModulus / (3. * (1. - 2. * fPoissonRatio));
    fShearModulus = fYoungModulus / (2. * (1. + fPoissonRatio));

    if (fDimension == 2) {
        fConstitutiveMatrix.resize(3,3);
        fConstitutiveMatrix.setZero();

        MatrixDouble I0 = MatrixDouble::Identity(3,3);
        I0(2,2) = 0.5;
        VecDouble m(3);
        m.setZero();
        m(0) = m(1) = 1.;

        fConstitutiveMatrix = 2. * fShearModulus * (I0 - 1./3. * m * m.transpose());

    } else if (fDimension == 3) {
        fConstitutiveMatrix.resize(6,6);
        fConstitutiveMatrix.setZero();

        MatrixDouble I0 = MatrixDouble::Identity(6,6);
        I0(3,3) = I0(4,4) = I0(5,5) = 0.5;
        VecDouble m(6);
        m.setZero();
        m(0) = m(1) = m(2) = 1.;

        fConstitutiveMatrix = 2. * fShearModulus * (I0 - 1./3. * m * m.transpose());
    } else {
        std::cerr << "Error: MixedElasticity only supports 2D and 3D problems." << std::endl;
        PanicButton();
    }
};

void MixedGlobalLocalEnrichment::ComputeStiffness(int &index, IntPointData &localdata, IntPointData &globaldata, MatrixDouble &Stiffness){
    
    //Only contribute stiffness if we are in a 2D element
    if (globaldata.fA0.rows() != 2) return;

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

void MixedGlobalLocalEnrichment::ComputeResidual(int &index, IntPointData &localdata, IntPointData &globaldata, VecDouble &Rhs){
   
    //Only contribute stiffness if we are in a 2D element
    if (globaldata.fA0.rows() == 2) return;

    int nphi = globaldata.fPhi.size();
    double WJ = localdata.fWeight * localdata.fJacA0;

    // MatrixDouble matB(3,2*nphi);
    // MatrixDouble matBEnr(3,2*nphi);
    // MatrixDouble matBTot(3, (matB.cols()+matBEnr.cols()));
    
    // matB.setZero();
    // matBEnr.setZero();
    // matBTot.setZero();

    auto force = fForceFunction;
    VecDouble forcingF(3);
    forcingF.setZero();
    VecDouble x_ = globaldata.fX;
    if (force) force(x_,forcingF);

    double uXInterp = localdata.fSol[0];
    double uYInterp = localdata.fSol[1];
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

    for (int i = nphi; i--; ){
        double shapeFi = globaldata.fPhi[i];
        //External force
        double Fx = forcingF[0] * uXInterp * shapeFi; 
        double Fy = forcingF[1] * uYInterp * shapeFi;
        Rhs[2*nphi + 2*i] += Fx * WJ;
        Rhs[2*nphi + 2*i+1] += Fy * WJ;
    };
};

int MixedGlobalLocalEnrichment::VariableIndex(const std::string &name) const{
    
    if(!strcmp("Displacement",name.c_str()))     return 1;
    
    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int MixedGlobalLocalEnrichment::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
        return 3;

    default:
        PanicButton();
        return -1;
    }
};

void MixedGlobalLocalEnrichment::Solution(IntPointData &data, int var, VecDouble &Sol) {

    //Solution
    if (var == 1){
        for (int i = 0; i < fNState; i++){
            Sol[i] = data.fSol[i];
        }
        return;
    };

}; 