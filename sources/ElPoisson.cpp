#include "ElPoisson.h"

void ElPoisson::ComputeStiffness(int &index, MatrixDouble &Stiffness){

    double WJ = fIntegData.fWeight * fIntegData.fJacA0 * getIntegPointWeightFunction(index);

    for (int i = Mesh()->NElNodes(); i-- ; ){       
        for (int j = Mesh()->NElNodes(); j-- ; ){            
            for (int k = Mesh()->Dimension(); k--;  ){
                //Diffusion matrix
                double K = fIntegData.fDPhiX0(i,k) * fIntegData.fDPhiX0(j,k);
                Stiffness(i,j) += K * WJ;
            }
        };
    };
    

}

void ElPoisson::ComputeResidual(int &index, VecDouble &Rhs){

    VecDouble fieldForce = Mesh()->getProblemParameters().GetFieldForce();
    auto force = Mesh()->getProblemParameters().getForcingFunction();
    int dim = Mesh()->Dimension();

    //Velocity Derivatives
    MatrixDouble du_dx(dim,dim);
    interpolateSolDerivatives(du_dx);

    double WJ = fIntegData.fWeight * fIntegData.fJacA0  * getIntegPointWeightFunction(index);

    VecDouble forcingF(1);
    VecDouble x_ = getIntegPointCoordinatesValue(index);
    if (force) force(x_,forcingF);

    for (int i = Mesh()->NElNodes(); i--; ){
        double shapeFi = Mesh()->getNumericalIntegration()-> phi_(i,index);

        //Viscosity
        double K = 0.;
        for (int l=dim; l--; ) K += fIntegData.fDPhiX0(i,l) * du_dx(0,l);

        //External force
        double F = (fieldForce[0] + forcingF[0]) * shapeFi;

        if(std::isnan(F)){
            PanicButton();
        }
        if(std::isnan(K)){
            PanicButton();
        }
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

    auto exactSol = Mesh()->getProblemParameters().getExactSolution();
    if (!exactSol) PanicButton();

    for(int it = 0; it < nQuad.getNumberOfIntegrationPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int i = DIM; i--; ) fIntegData.fAdimCoord[i] = nQuad.PointList(index,i);

        //Returns the quadrature integration weight
        fIntegData.fWeight = nQuad.WeightList(index);

        //Computes the jacobian matrix
        ComputeJacobian(index);
                    
        ComputeSpatialDerivatives();
        
        VecDouble uMEF_(1);
        interpolateSolution(index, uMEF_);
        MatrixDouble du_dxMEF(1,DIM);
        interpolateSolDerivatives(du_dxMEF);
        
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
        errors[0] += (u_[0]-uMEF_[0])*(u_[0]-uMEF_[0]) * fIntegData.fWeight * fIntegData.fJacA0 ;
        
        //Semi H1 state variable
        for (int m = DIM; m--; ){
            errors[1] += (gradU(m,0)-du_dxMEF(0,m))* (gradU(m,0)-du_dxMEF(0,m)) * fIntegData.fWeight * fIntegData.fJacA0;
        }

        index++;        
    }; 

    //H1 state variable
    errors[2] = errors[0]+errors[1];
}

void ElPoisson::ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs){

    for (int i = Mesh()->NElNodes(); i--; ){
        if ((Mesh()->NodeVec()[getConnectivity()[i]] -> getConstrains(0) == 1) ||
            (Mesh()->NodeVec()[getConnectivity()[i]] -> getConstrains(0) == 3))  {
            for (int j = Mesh()->NElNodes(); j--; ){
                Stiffness(i,j) = 0.;
                Stiffness(j,i) = 0.;
            };
            Stiffness(i,i) = 1.;
            Rhs[i] = 0.;
        }
    }
    // std::cout<<"Rhs -" << Rhs<<std::endl;

}