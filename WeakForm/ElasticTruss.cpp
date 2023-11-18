#include "ElasticTruss.h"

ElasticTruss::ElasticTruss(int matid, int dim, double young, double area) : WeakForm() {
    this->fMatId = matid;
    fDimension = dim;
    fNState = dim;
    fYoungModulus = young;
    fArea = area;
};


void ElasticTruss::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){

    data.fNeedsDSol = true;
    data.fDSolDx.resize(this->fDimension, this->fDimension);

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    int nphi = data.fPhi.size();
    double elementLenght = 2.*data.fJacA0;
    double K = fYoungModulus * fArea / elementLenght;



    MatrixDouble rotation(fDimension*nphi,fDimension*nphi);
    MatrixDouble matB(fDimension,fDimension*nphi);
    rotation.setZero();
    matB.setZero();
    double cosa = data.fAxes(0,0) / data.fJacA0;
    double sina = data.fAxes(1,0) / data.fJacA0;
    double check = sina*sina+cosa*cosa;
    for (int j = 0; j < nphi; j++){
        for (int i = 0; i < fDimension; i++){
            matB(i,fDimension*j + i) = data.fDPhiX0(j,0);
        }
        rotation(2*j  ,2*j  ) = cosa;
        rotation(2*j+1,2*j  ) = -sina;
        rotation(2*j  ,2*j+1) = sina;
        rotation(2*j+1,2*j+1) = cosa;
    }
    // std::cout << "rotation =\n"<< rotation << std::endl;
    Stiffness += rotation * matB.transpose() * matB * rotation.transpose() * WJ * K;

    // std::cout << "Stiffness =\n"<< Stiffness << std::endl;
}

void ElasticTruss::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    // int nphi = data.fPhi.size();

    // double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    // MatrixDouble matB(3,2*nphi);
    // matB.setZero();

    // auto force = fForceFunction;
    // VecDouble forcingF(fDimension);
    // VecDouble x_ = data.fX;
    // if (force) force(x_,forcingF);
    
    // for (int j = 0; j < nphi; j++){
    //     matB(0,fDimension*j  ) = data.fDPhiX0(j,0);
    //     matB(1,fDimension*j+1) = data.fDPhiX0(j,1);
    //     matB(2,fDimension*j  ) = data.fDPhiX0(j,1);
    //     matB(2,fDimension*j+1) = data.fDPhiX0(j,0);
    // }
    
    // VecDouble strain(3);
    // strain.setZero();
    // strain[0] = data.fDSolDx(0,0);
    // strain[1] = data.fDSolDx(1,1);
    // strain[2] = data.fDSolDx(0,1)+data.fDSolDx(1,0);

    // Rhs -= matB.transpose() * fConstitutiveMatrix * strain * WJ;

    // for (int i = nphi; i--; ){
    //     double shapeFi = data.fPhi[i];
    //     //External force
    //     double Fx = forcingF[0] * shapeFi;
    //     double Fy = forcingF[1] * shapeFi;
    //     Rhs[2*i  ] += Fx * WJ;
    //     Rhs[2*i+1] += Fy * WJ;
    // };

    
};

void ElasticTruss::ComputeError(IntPointData &data, VecDouble &errors){
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

int ElasticTruss::VariableIndex(const std::string &name) const{
    
    if(!strcmp("Displacement",name.c_str()))           return 1;
    if(!strcmp("NormalStress",name.c_str()))           return 2;
    if(!strcmp("ExactDisplacement",name.c_str()))      return 3;
    if(!strcmp("ExactNormalStress",name.c_str()))      return 4;
    if(!strcmp("ExactForce",name.c_str()))             return 5;

    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int ElasticTruss::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
    case 3:
    case 5:
        return 3;
    case 2:
    case 4:
        return 1;

    default:
        PanicButton();
        return -1;
    }
};

void ElasticTruss::Solution(IntPointData &data, int var, VecDouble &Sol) {

    //Displacement
    if (var == 1){
        Sol[0] = data.fSol[0];
        Sol[1] = data.fSol[1];
        Sol[2] = 0.;
        return;
    };

    //NormalStress
    if (var == 2){
        double epsilon = data.fDSolDx.norm();
        Sol[0] = fYoungModulus * epsilon;
        return;
    };

    VecDouble forcingF(fDimension);
    VecDouble x_ = data.fX;
    if (fForceFunction) fForceFunction(x_,forcingF);

    VecDouble disp(fDimension);
    MatrixDouble gradDisp(fDimension,3);
    if (fExactSol) fExactSol(x_,disp,gradDisp);

    //Exact Displacement
    if (var == 3){
        Sol[0] = disp[0];
        Sol[1] = disp[1];
        Sol[2] = 0.;
        return;
    };
    
    //Exact Sigma X
    if (var == 4){
        double epsilon = gradDisp.norm();
        Sol[0] = fYoungModulus * epsilon;
        return;
    };

    //Exact Force
    if (var == 5){
        Sol[0] = forcingF[0];
        Sol[1] = forcingF[1];
        Sol[2] = 0.;
        return;
    };

}; 