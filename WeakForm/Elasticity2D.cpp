#include "Elasticity2D.h"

void Elasticity2D::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){

    data.fNeedsDSol = true;
    data.fDSolDx.resize(this->fDimension, this->fDimension);

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    int nphi = data.fPhi.size();
    MatrixDouble matD(3,2*nphi);
    matD.setZero();

    MatrixDouble Hooke(3,3);
    Hooke.setZero();
    // For EPT
    if (fPlaneStress){
        double k = fYoungModulus / (1. - fPoissonRatio * fPoissonRatio);
        Hooke(0,0) = k;
        Hooke(0,1) = k * fPoissonRatio;
        Hooke(1,0) = k * fPoissonRatio;
        Hooke(1,1) = k;
        Hooke(2,2) = k * (1. - fPoissonRatio) * 0.5;
    } else {
        std::cout << "Please Implement me" <<std::endl;
    }
    

    for (int j = 0; j < nphi; j++){
        matD(0,fDimension*j  ) = data.fDPhiX0(j,0);
        matD(1,fDimension*j+1) = data.fDPhiX0(j,1);
        matD(2,fDimension*j  ) = data.fDPhiX0(j,1);
        matD(2,fDimension*j+1) = data.fDPhiX0(j,0);
    }
    
    Stiffness += matD.transpose() * Hooke * matD * WJ;

    // std::cout << "Stiffness =\n"<< Stiffness << std::endl;
}

void Elasticity2D::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    int nphi = data.fPhi.size();

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    MatrixDouble matD(3,2*nphi);
    matD.setZero();

    auto force = fForceFunction;
    VecDouble forcingF(fDimension);
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    MatrixDouble Hooke(3,3);
    Hooke.setZero();
    // For EPT
    if (fPlaneStress){
        double k = fYoungModulus / (1. - fPoissonRatio * fPoissonRatio);
        Hooke(0,0) = k;
        Hooke(0,1) = k * fPoissonRatio;
        Hooke(1,0) = k * fPoissonRatio;
        Hooke(1,1) = k;
        Hooke(2,2) = k * (1. - fPoissonRatio) * 0.5;
    } else {
        std::cout << "Please Implement me" <<std::endl;
    }
    
    for (int j = 0; j < nphi; j++){
        matD(0,fDimension*j  ) = data.fDPhiX0(j,0);
        matD(1,fDimension*j+1) = data.fDPhiX0(j,1);
        matD(2,fDimension*j  ) = data.fDPhiX0(j,1);
        matD(2,fDimension*j+1) = data.fDPhiX0(j,0);
    }
    
    VecDouble strain(3);
    strain.setZero();
    strain[0] = data.fDSolDx(0,0);
    strain[1] = data.fDSolDx(1,1);
    strain[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);

    Rhs -= matD.transpose() * Hooke * strain * WJ;

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

