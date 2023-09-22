#include "ElPoisson.h"

void ElPoisson::ComputeStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, MatrixDouble &Stiffness){

    double WJ = weight_ * djac_ * getIntegPointWeightFunction(index);

    for (int i = Mesh()->nElNodes; i-- ; ){       
        for (int j = Mesh()->nElNodes; j-- ; ){            
            for (int k = Mesh()->Dimension(); k--;  ){
                //Diffusion matrix
                double K = dphi_dx(i,k) * dphi_dx(j,k);
                Stiffness(i,j) += K * WJ;
            }
        };
    };


}

void ElPoisson::ComputeResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, VecDouble &Rhs){

    VecDouble fieldForce = Mesh()->getProblemParameters().GetFieldForce();
    auto force = Mesh()->getProblemParameters().getForcingFunction();
    int dim = Mesh()->Dimension();

    //Velocity Derivatives
    MatrixDouble du_dx(dim,dim), duprev_dx(dim,dim);
    interpolateVelDerivatives(dphi_dx, du_dx, duprev_dx);

    double WJ = weight_ * djac_  * getIntegPointWeightFunction(index);

    VecDouble forcingF(1);
    VecDouble x_ = getIntegPointCoordinatesValue(index);
    if (force) force(x_,forcingF);

    for (int i = Mesh()->nElNodes; i--; ){
        double shapeFi = Mesh()->getNumericalIntegration()-> phi_(i,index);

        //Viscosity
        double K = 0.;
        for (int l=dim; l--; ) K += dphi_dx(i,l) * du_dx(0,l);

        //External force
        double F = (fieldForce[0] + forcingF[0]) * shapeFi;
        Rhs[i] += (-K + F) * WJ;
    };
};

void ElPoisson::ComputeError(VecDouble &errors){
    int index = 0;
    errors.resize(3);
    errors.setZero();
    int DIM = Mesh()->Dimension();
    int DEG = Mesh()->GetDefaultOrder();

    IntegQuadrature nQuad(DIM,DEG);
    ShapeFunction shapeQuad(DIM,DEG);

    MatrixDouble dphi_dx(Mesh()->nElNodes,DIM);
    MatrixDouble ainv_(DIM,DIM);
    VecDouble xsi(DIM);
    double weight_;

    auto exactSol = Mesh()->getProblemParameters().getExactSolution();
    if (!exactSol) PanicButton();

    for(int it = 0; it < nQuad.getNumberOfIntegrationPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int i = DIM; i--; ) xsi[i] = nQuad.PointList(index,i);

        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        double djac_ = 0.;
        //Computes the jacobian matrix
        getJacobianMatrix(xsi, ainv_, djac_, index);
                    
        getSpatialDerivatives(xsi, ainv_, dphi_dx);
        
        VecDouble uMEF_(DIM), uPrev_(DIM);
        interpolateVelocity(index, uMEF_, uPrev_);
        MatrixDouble du_dxMEF(DIM,DIM), duprev_dx(DIM,DIM);
        interpolateVelDerivatives(dphi_dx, du_dxMEF, duprev_dx);
        
        VecDouble u_(1);
        MatrixDouble gradU(DIM,1);

        VecDouble xna_ = getIntegPointCoordinatesValue(index);
        
        exactSol(xna_,u_,gradU);


        //Consider Arlequin weight function
        u_ *= getIntegPointWeightFunction(index);
        gradU *= getIntegPointWeightFunction(index);
        uMEF_ *= getIntegPointWeightFunction(index);
        du_dxMEF *= getIntegPointWeightFunction(index);

        //L2 state variable
        errors[0] += (u_[0]-uMEF_[0])*(u_[0]-uMEF_[0]) * weight_ * djac_ ;
        
        //Semi H1 state variable
        for (int m = DIM; m--; ){
            errors[1] += (gradU(m,0)-du_dxMEF(0,m))* (gradU(m,0)-du_dxMEF(0,m)) * weight_ * djac_;
        }

        index++;        
    }; 

    //H1 state variable
    errors[2] = errors[0]+errors[1];
}

void ElPoisson::ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs){

    for (int i = Mesh()->nElNodes; i--; ){
        if ((Mesh()->NodeVec()[getConnectivity()[i]] -> getConstrains(0) == 1) ||
            (Mesh()->NodeVec()[getConnectivity()[i]] -> getConstrains(0) == 3))  {
            for (int j = Mesh()->nElNodes; j--; ){
                Stiffness(i,j) = 0.;
                Stiffness(j,i) = 0.;
            };
            Stiffness(i,i) = 1.;
            Rhs[i] = Mesh()->NodeVec()[getConnectivity()[i]]->GetSolution(0);
        }
    }
    // std::cout<<"Rhs -" << Rhs<<std::endl;

}