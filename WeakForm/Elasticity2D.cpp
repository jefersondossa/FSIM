#include "Elasticity2D.h"

Elasticity2D::Elasticity2D(int matid, double young, double poisson, bool planes) : WeakForm() {
    this->fMatId = matid;
    fDimension = 2;
    fNState = 2;
    fYoungModulus = young;
    fPoissonRatio = poisson;
    fPlaneStress = planes;
    fConstitutiveMatrix.resize(3,3);
    fConstitutiveMatrix.setZero();
    if (fPlaneStress){//Plane Stress Matrix
        double k = fYoungModulus / (1. - fPoissonRatio * fPoissonRatio);
        fConstitutiveMatrix(0,0) = k;
        fConstitutiveMatrix(0,1) = k * fPoissonRatio;
        fConstitutiveMatrix(1,0) = k * fPoissonRatio;
        fConstitutiveMatrix(1,1) = k;
        fConstitutiveMatrix(2,2) = k * (1. - fPoissonRatio) * 0.5;
    } else {//Plane Strain Matrix
        double G = fYoungModulus / (2. * ( 1. + fPoissonRatio));
        double k = 2.*G / (1.-2.*fPoissonRatio);
        fConstitutiveMatrix(0,0) = (1.-fPoissonRatio) * k;
        fConstitutiveMatrix(0,1) = k * fPoissonRatio;
        fConstitutiveMatrix(1,0) = k * fPoissonRatio;
        fConstitutiveMatrix(1,1) = (1.-fPoissonRatio) * k;
        fConstitutiveMatrix(2,2) = 2. * G;
    }
};


void Elasticity2D::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){

    if (!data.fNeedsDSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
    }

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    int nphi = data.fPhi.size();
    MatrixDouble matB(3,2*nphi);
    matB.setZero();

    for (int j = 0; j < nphi; j++){
        matB(0,2*j  ) = data.fDPhiX0(j,0);
        matB(1,2*j+1) = data.fDPhiX0(j,1);
        matB(2,2*j  ) = data.fDPhiX0(j,1);
        matB(2,2*j+1) = data.fDPhiX0(j,0);
    }
    
    Stiffness += matB.transpose() * fConstitutiveMatrix * matB * WJ;

    // std::cout << "Stiffness =\n"<< Stiffness << std::endl;
}

void Elasticity2D::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    int nphi = data.fPhi.size();

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    MatrixDouble matB(3,2*nphi);
    matB.setZero();

    auto force = fForceFunction;
    VecDouble forcingF(fDimension);
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);
    
    for (int j = 0; j < nphi; j++){
        matB(0,fDimension*j  ) = data.fDPhiX0(j,0);
        matB(1,fDimension*j+1) = data.fDPhiX0(j,1);
        matB(2,fDimension*j  ) = data.fDPhiX0(j,1);
        matB(2,fDimension*j+1) = data.fDPhiX0(j,0);
    }
    
    VecDouble strain(3);
    strain.setZero();
    strain[0] = data.fDSolDx(0,0);
    strain[1] = data.fDSolDx(1,1);
    strain[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);

    Rhs -= matB.transpose() * fConstitutiveMatrix * strain * WJ;

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
    // int index = 0;
    // errors.resize(3);
    // errors.setZero();
    // int DIM = this->Mesh()->Dimension();
    // int DEG = this->Mesh()->GetDefaultOrder();

    // IntegQuadrature nQuad(DIM,DEG);

    // MatrixDouble dphi_dx(tshape::NElNodes,DIM);
    // MatrixDouble ainv_(DIM,DIM);
    // VecDouble xsi(DIM);

    // auto exactSol = this->Mesh()->getProblemParameters().getExactSolution();
    // if (!exactSol) PanicButton();

    // for(int it = 0; it < nQuad.getNumberOfIntegrationPoints(); it++){

    //     //Defines the integration points adimentional coordinates
    //     for (int i = DIM; i--; ) xsi[i] = nQuad.PointList(index,i);

    //     //Returns the quadrature integration weight
    //     data.fWeight = nQuad.WeightList(index);

    //     //Computes the jacobian matrix
    //     this->ComputeJacobian(index);
                    
    //     this->ComputeSpatialDerivatives();
        
    //     VecDouble uMEF_(DIM);
    //     this->interpolateSolution(index, uMEF_);
    //     MatrixDouble du_dxMEF(DIM,DIM);
    //     this->interpolateSolDerivatives(du_dxMEF);
        
    //     VecDouble u_(DIM);
    //     MatrixDouble gradU(DIM,DIM);

    //     VecDouble x_ = this->getIntegPointCoordinatesValue(index);
        
    //     exactSol(x_,u_,gradU);

    //     //Consider Arlequin weight function
    //     u_ *= this->getIntegPointWeightFunction(index);
    //     gradU *= this->getIntegPointWeightFunction(index);
    //     uMEF_ *= this->getIntegPointWeightFunction(index);
    //     du_dxMEF *= this->getIntegPointWeightFunction(index);

    //     //L2 displacement
    //     errors[0] += ((u_[0]-uMEF_[0])*(u_[0]-uMEF_[0]) + (u_[1]-uMEF_[1])*(u_[1]-uMEF_[1]))
    //                   * data.fWeight * data.fJacA0 ;
        

    //     // std::cout << "Stress and Energy norms not implemented yet\n";
    //     // //Semi H1 state variable
    //     // for (int m = DIM; m--; ){
    //     //     errors[1] += (gradU[m]-du_dxMEF(0,m))* (gradU[m]-du_dxMEF(0,m)) * weight_ * djac_;
    //     // }

    //     index++;        
    // }; 

    // //H1 state variable
    // errors[2] = errors[0]+errors[1];
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

    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int Elasticity2D::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
    case 8:
    case 15:
    case 16:
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

}; 


MatrixDouble Elasticity2D::ConstitutiveMatrix(){
    return fConstitutiveMatrix;
}