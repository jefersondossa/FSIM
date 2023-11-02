#include "ElElasticity2D.h"

template<class tshape>
void ElElasticity2D<tshape>::ComputeStiffness(int &index, MatrixDouble &Stiffness){

    double WJ = this->fIntegData.fWeight * this->fIntegData.fJacA0 * this->getIntegPointWeightFunction(index);
    MatrixDouble matD(3,2*tshape::NElNodes);
    matD.setZero();
    int DIM = this->Mesh()->Dimension();

    MatrixDouble Hooke(3,3);
    Hooke.setZero();
    // For EPT
    double elastic_ = this->Mesh()->getProblemParameters().GetYoungModulus();
    double poisson_ = this->Mesh()->getProblemParameters().GetPoissonRatio();
    double k = elastic_ / (1. - poisson_ * poisson_);
    Hooke(0,0) = k;
    Hooke(0,1) = k * poisson_;
    Hooke(1,0) = k * poisson_;
    Hooke(1,1) = k;
    Hooke(2,2) = k * (1. - poisson_) * 0.5;

    for (int j = 0; j < tshape::NElNodes; j++){
        matD(0,DIM*j  ) = this->fIntegData.fDPhiX0(j,0);
        matD(1,DIM*j+1) = this->fIntegData.fDPhiX0(j,1);
        matD(2,DIM*j  ) = this->fIntegData.fDPhiX0(j,1);
        matD(2,DIM*j+1) = this->fIntegData.fDPhiX0(j,0);
    }
    
    Stiffness += matD.transpose() * Hooke * matD * WJ;


    // std::cout << "Stiffness =\n"<< Stiffness << std::endl;
}

template<class tshape>
void ElElasticity2D<tshape>::ComputeResidual(int &index, VecDouble &Rhs){

    VecDouble fieldForce = this->Mesh()->getProblemParameters().GetFieldForce();
    auto force = this->Mesh()->getProblemParameters().getForcingFunction();
    int DIM = this->Mesh()->Dimension();

    //Velocity Derivatives
    MatrixDouble du_dx(DIM,DIM);
    this->interpolateSolDerivatives(du_dx);

    double WJ = this->fIntegData.fWeight * this->fIntegData.fJacA0 * this->getIntegPointWeightFunction(index);
    MatrixDouble matD(3,2*tshape::NElNodes);
    matD.setZero();

    MatrixDouble Hooke(3,3);
    Hooke.setZero();
    // For EPT
    double elastic_ = this->Mesh()->getProblemParameters().GetYoungModulus();
    double poisson_ = this->Mesh()->getProblemParameters().GetPoissonRatio();
    double k = elastic_ / (1. - poisson_ * poisson_);
    Hooke(0,0) = k;
    Hooke(0,1) = k * poisson_;
    Hooke(1,0) = k * poisson_;
    Hooke(1,1) = k;
    Hooke(2,2) = k * (1. - poisson_) * 0.5;

    for (int j = 0; j < tshape::NElNodes; j++){
        matD(0,DIM*j  ) = this->fIntegData.fDPhiX0(j,0);
        matD(1,DIM*j+1) = this->fIntegData.fDPhiX0(j,1);
        matD(2,DIM*j  ) = this->fIntegData.fDPhiX0(j,1);
        matD(2,DIM*j+1) = this->fIntegData.fDPhiX0(j,0);
    }
    
    VecDouble strain(3);
    strain.setZero();
    strain[0] = du_dx(0,0);
    strain[1] = du_dx(1,1);
    strain[2] = du_dx(0,1)+du_dx(1,0);

    Rhs -= matD.transpose() * Hooke * strain * WJ;

    VecDouble forcingF(2);
    VecDouble x_ = this->getIntegPointCoordinatesValue(index);
    if (force) force(x_,forcingF);

    for (int i = tshape::NElNodes; i--; ){
        double shapeFi = this->fIntegData.fPhi[i];
        //External force
        double Fx = (fieldForce[0] + forcingF[0]) * shapeFi;
        double Fy = (fieldForce[1] + forcingF[1]) * shapeFi;
        Rhs[2*i  ] += Fx * WJ;
        Rhs[2*i+1] += Fy * WJ;
    };

    
};
template<class tshape>
void ElElasticity2D<tshape>::ComputeError(VecDouble &errors){
    int index = 0;
    errors.resize(3);
    errors.setZero();
    int DIM = this->Mesh()->Dimension();
    int DEG = this->Mesh()->GetDefaultOrder();

    IntegQuadrature nQuad(DIM,DEG);

    MatrixDouble dphi_dx(tshape::NElNodes,DIM);
    MatrixDouble ainv_(DIM,DIM);
    VecDouble xsi(DIM);

    auto exactSol = this->Mesh()->getProblemParameters().getExactSolution();
    if (!exactSol) PanicButton();

    for(int it = 0; it < nQuad.getNumberOfIntegrationPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int i = DIM; i--; ) xsi[i] = nQuad.PointList(index,i);

        //Returns the quadrature integration weight
        this->fIntegData.fWeight = nQuad.WeightList(index);

        //Computes the jacobian matrix
        this->ComputeJacobian(index);
                    
        this->ComputeSpatialDerivatives();
        
        VecDouble uMEF_(DIM);
        this->interpolateSolution(index, uMEF_);
        MatrixDouble du_dxMEF(DIM,DIM);
        this->interpolateSolDerivatives(du_dxMEF);
        
        VecDouble u_(DIM);
        MatrixDouble gradU(DIM,DIM);

        VecDouble x_ = this->getIntegPointCoordinatesValue(index);
        
        exactSol(x_,u_,gradU);

        //Consider Arlequin weight function
        u_ *= this->getIntegPointWeightFunction(index);
        gradU *= this->getIntegPointWeightFunction(index);
        uMEF_ *= this->getIntegPointWeightFunction(index);
        du_dxMEF *= this->getIntegPointWeightFunction(index);

        //L2 displacement
        errors[0] += ((u_[0]-uMEF_[0])*(u_[0]-uMEF_[0]) + (u_[1]-uMEF_[1])*(u_[1]-uMEF_[1]))
                      * this->fIntegData.fWeight * this->fIntegData.fJacA0 ;
        

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

template<class tshape>
void ElElasticity2D<tshape>::ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs){

    for (int i = tshape::NElNodes; i--; ){
        int nstate = this->Mesh()->NodeVec()[this->getConnectivity()[i]]->GetNStateVariables();
        for (int istate = 0; istate < nstate; istate++){
            if ((this->Mesh()->NodeVec()[this->getConnectivity()[i]] -> getConstrains(istate) == 1) ||
                (this->Mesh()->NodeVec()[this->getConnectivity()[i]] -> getConstrains(istate) == 3))  {
                for (int j = tshape::NElNodes*nstate; j--; ){
                    Stiffness(nstate*i+istate,j) = 0.;
                    Stiffness(j,nstate*i+istate) = 0.;
                };
                Stiffness(nstate*i+istate,nstate*i+istate) = 1.;
                Rhs[nstate*i+istate] = 0.;
            }
        }
    }

}



#include "ShapeHexahedron.h"
#include "ShapeOneD.h"
#include "ShapeQuadrilateralLin.h"
#include "ShapePoint.h"
#include "ShapeTetrahedronLin.h"
#include "ShapeTetrahedronQua.h"
#include "ShapeTetrahedronCub.h"
#include "ShapeTriangleLin.h"
#include "ShapeTriangleQua.h"
#include "ShapeTriangleCub.h"

template class ElElasticity2D<ShapePoint>;
template class ElElasticity2D<ShapeOneD>;
template class ElElasticity2D<ShapeTriangleLin>;
template class ElElasticity2D<ShapeTriangleQua>;
template class ElElasticity2D<ShapeTriangleCub>;
template class ElElasticity2D<ShapeQuadrilateralLin>;
template class ElElasticity2D<ShapeTetrahedronLin>;
template class ElElasticity2D<ShapeTetrahedronQua>;
template class ElElasticity2D<ShapeTetrahedronCub>;
template class ElElasticity2D<ShapeHexahedron>;