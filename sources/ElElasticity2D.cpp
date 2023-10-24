#include "ElElasticity2D.h"

void ElElasticity2D::ComputeStiffness(int &index, MatrixDouble &Stiffness){

    double WJ = fIntegData.fWeight * fIntegData.fJacA0 * getIntegPointWeightFunction(index);
    MatrixDouble matD(3,2*Mesh()->NElNodes());
    matD.setZero();
    int DIM = Mesh()->Dimension();

    MatrixDouble Hooke(3,3);
    Hooke.setZero();
    // For EPT
    double elastic_ = Mesh()->getProblemParameters().GetYoungModulus();
    double poisson_ = Mesh()->getProblemParameters().GetPoissonRatio();
    double k = elastic_ / (1. - poisson_ * poisson_);
    Hooke(0,0) = k;
    Hooke(0,1) = k * poisson_;
    Hooke(1,0) = k * poisson_;
    Hooke(1,1) = k;
    Hooke(2,2) = k * (1. - poisson_) * 0.5;

    for (int j = 0; j < Mesh()->NElNodes(); j++){
        matD(0,DIM*j  ) = fIntegData.fDPhiX0(j,0);
        matD(1,DIM*j+1) = fIntegData.fDPhiX0(j,1);
        matD(2,DIM*j  ) = fIntegData.fDPhiX0(j,1);
        matD(2,DIM*j+1) = fIntegData.fDPhiX0(j,0);
    }
    
    Stiffness += matD.transpose() * Hooke * matD * WJ;


    // std::cout << "Stiffness =\n"<< Stiffness << std::endl;
}

void ElElasticity2D::ComputeResidual(int &index, VecDouble &Rhs){

    VecDouble fieldForce = Mesh()->getProblemParameters().GetFieldForce();
    auto force = Mesh()->getProblemParameters().getForcingFunction();
    int DIM = Mesh()->Dimension();

    //Velocity Derivatives
    MatrixDouble du_dx(DIM,DIM);
    interpolateSolDerivatives(du_dx);

    double WJ = fIntegData.fWeight * fIntegData.fJacA0 * getIntegPointWeightFunction(index);
    MatrixDouble matD(3,2*Mesh()->NElNodes());
    matD.setZero();

    MatrixDouble Hooke(3,3);
    Hooke.setZero();
    // For EPT
    double elastic_ = Mesh()->getProblemParameters().GetYoungModulus();
    double poisson_ = Mesh()->getProblemParameters().GetPoissonRatio();
    double k = elastic_ / (1. - poisson_ * poisson_);
    Hooke(0,0) = k;
    Hooke(0,1) = k * poisson_;
    Hooke(1,0) = k * poisson_;
    Hooke(1,1) = k;
    Hooke(2,2) = k * (1. - poisson_) * 0.5;

    for (int j = 0; j < Mesh()->NElNodes(); j++){
        matD(0,DIM*j  ) = fIntegData.fDPhiX0(j,0);
        matD(1,DIM*j+1) = fIntegData.fDPhiX0(j,1);
        matD(2,DIM*j  ) = fIntegData.fDPhiX0(j,1);
        matD(2,DIM*j+1) = fIntegData.fDPhiX0(j,0);
    }
    
    VecDouble strain(3);
    strain.setZero();
    strain[0] = du_dx(0,0);
    strain[1] = du_dx(1,1);
    strain[2] = du_dx(0,1)+du_dx(1,0);

    Rhs -= matD.transpose() * Hooke * strain * WJ;

    VecDouble forcingF(2);
    VecDouble x_ = getIntegPointCoordinatesValue(index);
    if (force) force(x_,forcingF);

    for (int i = Mesh()->NElNodes(); i--; ){
        double shapeFi = Mesh()->getNumericalIntegration()-> phi_(i,index);
        //External force
        double Fx = (fieldForce[0] + forcingF[0]) * shapeFi;
        double Fy = (fieldForce[1] + forcingF[1]) * shapeFi;
        Rhs[2*i  ] += Fx * WJ;
        Rhs[2*i+1] += Fy * WJ;
    };

    
};

void ElElasticity2D::ComputeError(VecDouble &errors){
    int index = 0;
    errors.resize(3);
    errors.setZero();
    int DIM = Mesh()->Dimension();
    int DEG = Mesh()->GetDefaultOrder();

    IntegQuadrature nQuad(DIM,DEG);
    ShapeFunction shapeQuad(DIM,DEG);

    MatrixDouble dphi_dx(Mesh()->NElNodes(),DIM);
    MatrixDouble ainv_(DIM,DIM);
    VecDouble xsi(DIM);

    auto exactSol = Mesh()->getProblemParameters().getExactSolution();
    if (!exactSol) PanicButton();

    for(int it = 0; it < nQuad.getNumberOfIntegrationPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int i = DIM; i--; ) xsi[i] = nQuad.PointList(index,i);

        //Returns the quadrature integration weight
        fIntegData.fWeight = nQuad.WeightList(index);

        //Computes the jacobian matrix
        ComputeJacobian(index);
                    
        ComputeSpatialDerivatives();
        
        VecDouble uMEF_(DIM);
        interpolateSolution(index, uMEF_);
        MatrixDouble du_dxMEF(DIM,DIM);
        interpolateSolDerivatives(du_dxMEF);
        
        VecDouble u_(DIM);
        MatrixDouble gradU(DIM,DIM);

        VecDouble x_ = getIntegPointCoordinatesValue(index);
        
        exactSol(x_,u_,gradU);

        //Consider Arlequin weight function
        u_ *= getIntegPointWeightFunction(index);
        gradU *= getIntegPointWeightFunction(index);
        uMEF_ *= getIntegPointWeightFunction(index);
        du_dxMEF *= getIntegPointWeightFunction(index);

        //L2 displacement
        errors[0] += ((u_[0]-uMEF_[0])*(u_[0]-uMEF_[0]) + (u_[1]-uMEF_[1])*(u_[1]-uMEF_[1]))
                      * fIntegData.fWeight * fIntegData.fJacA0 ;
        

        // std::cout << "Stress and Energy norms not implemented yet\n";
        // //Semi H1 state variable
        // for (int m = DIM; m--; ){
        //     errors[1] += (gradU[m]-du_dxMEF(0,m))* (gradU[m]-du_dxMEF(0,m)) * weight_ * djac_;
        // }

        index++;        
    }; 

    //H1 state variable
    errors[2] = errors[0]+errors[1];
}

void ElElasticity2D::ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs){

    for (int i = Mesh()->NElNodes(); i--; ){
        int nstate = Mesh()->NodeVec()[getConnectivity()[i]]->GetNStateVariables();
        for (int istate = 0; istate < nstate; istate++){
            if ((Mesh()->NodeVec()[getConnectivity()[i]] -> getConstrains(istate) == 1) ||
                (Mesh()->NodeVec()[getConnectivity()[i]] -> getConstrains(istate) == 3))  {
                for (int j = Mesh()->NElNodes()*nstate; j--; ){
                    Stiffness(nstate*i+istate,j) = 0.;
                    Stiffness(j,nstate*i+istate) = 0.;
                };
                Stiffness(nstate*i+istate,nstate*i+istate) = 1.;
                Rhs[nstate*i+istate] = 0.;
            }
        }
    }

}
