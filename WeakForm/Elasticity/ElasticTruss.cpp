#include "ElasticTruss.h"
#include "Element.h"

ElasticTruss::ElasticTruss(int matid, int dim, double young, double area) : WeakForm() {
    this->fMatId = matid;
    fDimension = dim;
    fNState = dim;
    fYoungModulus = young;
    fArea = area;
    this->fType = WeakFormType::kElasticTruss;
};


void ElasticTruss::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    if (!data.fNeedsDSol || !data.fNeedsSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
        data.fNeedsSol = true;
        data.fSol.resize(fNState);
    }
    
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    int nphi = data.fPhi.size();
    double elementLenght = 2.*data.fJacA0;
    double K = fYoungModulus * fArea / elementLenght;


    MatrixDouble rotation(fDimension*nphi,fDimension*nphi);
    MatrixDouble matB(fDimension,fDimension*nphi);
    rotation.setZero();
    matB.setZero();
    double cosa = data.fAxes0(0,0) / data.fJacA0;
    double sina = data.fAxes0(1,0) / data.fJacA0;
    double check = sina*sina+cosa*cosa;
    for (int j = 0; j < nphi; j++){
        // for (int i = 0; i < fDimension; i++){
            matB(0,fDimension*j) = data.fDPhiX0(0,j);
        // }
        rotation(2*j  ,2*j  ) = cosa;
        rotation(2*j+1,2*j  ) = sina;
        rotation(2*j  ,2*j+1) = -sina;
        rotation(2*j+1,2*j+1) = cosa;
    }
    // std::cout << "matB =\n"<< matB << std::endl;
    // std::cout << "rotation =\n"<< rotation << std::endl;
    // std::cout << "K =\n"<< rotation * matB.transpose() * matB * rotation.transpose() << std::endl;
    Stiffness += rotation * matB.transpose() * matB * rotation.transpose() * WJ * elementLenght * K;

}

void ElasticTruss::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    int nphi = data.fPhi.size();
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    auto force = fForceFunction;
    VecDouble forcingF(fDimension);
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);
    
    double elementLenght = 2.*data.fJacA0;
    double K = fYoungModulus * fArea / elementLenght;

    MatrixDouble rotation(fDimension*nphi,fDimension*nphi);
    MatrixDouble matB(fDimension,fDimension*nphi);
    rotation.setZero();
    matB.setZero();
    double cosa = data.fAxes0(0,0) / data.fJacA0;
    double sina = data.fAxes0(1,0) / data.fJacA0;
    double check = sina*sina+cosa*cosa;
    for (int j = 0; j < nphi; j++){
        // for (int i = 0; i < fDimension; i++){
            matB(0,fDimension*j) = data.fDPhiX0(0,j);
        // }
        rotation(2*j  ,2*j  ) = cosa;
        rotation(2*j+1,2*j  ) = sina;
        rotation(2*j  ,2*j+1) = -sina;
        rotation(2*j+1,2*j+1) = cosa;
    }

    VecDouble sol(2);// = data.fSol;
    sol[0] = -data.fDSolDx(1,0)*sina;
    sol[1] = +data.fDSolDx(0,0)*cosa;

    // std::cout << "rotation =\n"<< rotation << std::endl;
    Rhs += rotation * matB.transpose() * sol * WJ * elementLenght * K;

    for (int i = nphi; i--; ){
        double shapeFi = data.fPhi[i];
        //External force
        double Fx = forcingF[0] * shapeFi;
        double Fy = forcingF[1] * shapeFi;
        Rhs[2*i  ] += Fx * WJ;
        Rhs[2*i+1] += Fy * WJ;
    };
    // std::cout << "RHS = " << Rhs << std::endl;
    
};

void ElasticTruss::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();

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
    if(!strcmp("Stress",name.c_str()))           return 2;
    if(!strcmp("Strain",name.c_str()))           return 6;
    if(!strcmp("AxialForce",name.c_str()))           return 7;
    if(!strcmp("ExactDisplacement",name.c_str()))      return 3;
    
    // std::cout << "Post Process variable not implemented \n";
    // PanicButton();
    return -1;
};

int ElasticTruss::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
    case 3:
        return 3;
    case 2:
    case 6:
    case 7:
        return 1;

    default:
        // PanicButton();
        return -1;
    }
};

void ElasticTruss::Solution(IntPointData &data, int var, VecDouble &Sol) {

    double cosa = data.fAxes0(0,0) / data.fJacA0;
    double sina = data.fAxes0(1,0) / data.fJacA0;

    VecDouble nodalSol(4);
    auto connect = data.fElement->getConnectivity();
    for (int i = 0; i < 2; i++){
        nodalSol[2*i  ] = connect[i]->GetSolution(0);
        nodalSol[2*i+1] = connect[i]->GetSolution(1);
    }
    double elementLenght = 2.*data.fJacA0;

    double strain = (nodalSol[2]-nodalSol[0])*cosa/elementLenght + (nodalSol[3]-nodalSol[1])*sina/elementLenght;

    //Displacement
    if (var == 1){
        Sol[0] = data.fSol[0];
        Sol[1] = data.fSol[1];
        Sol[2] = 0.;
        return;
    };

    //NormalStress
    if (var == 2){
        Sol[0] = fYoungModulus * strain;
        return;
    };
    // Strain
    if (var == 6){
        Sol[0] = strain;
        return;
    }

    // Axial Force
    if (var == 7){
        Sol[0] = fYoungModulus * fArea * strain;
        return;
    }

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
}; 

MatrixDouble &ElasticTruss::ConstitutiveMatrix(){
    MatrixDouble constitutive(1,1);
    constitutive(0,0) = fYoungModulus;
    return constitutive;
}