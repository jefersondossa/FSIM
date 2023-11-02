#include "ElPoisson.h"

template <class tshape>
void ElPoisson<tshape>::ComputeStiffness(int &index, MatrixDouble &Stiffness){

    double WJ = this->fIntegData.fWeight * this->fIntegData.fJacA0 * this->getIntegPointWeightFunction(index);

    for (int i = tshape::NElNodes; i-- ; ){       
        for (int j = tshape::NElNodes; j-- ; ){            
            for (int k = this->Mesh()->Dimension(); k--;  ){
                //Diffusion matrix
                double K = this->fIntegData.fDPhiX0(i,k) * this->fIntegData.fDPhiX0(j,k);
                Stiffness(i,j) += K * WJ;
            }
        };
    };
}

template <class tshape>
void ElPoisson<tshape>::ComputeResidual(int &index, VecDouble &Rhs){

    VecDouble fieldForce = this->Mesh()->getProblemParameters().GetFieldForce();
    auto force = this->Mesh()->getProblemParameters().getForcingFunction();
    int dim = this->Mesh()->Dimension();

    //Velocity Derivatives
    MatrixDouble du_dx(dim,dim);
    this->interpolateSolDerivatives(du_dx);

    double WJ = this->fIntegData.fWeight * this->fIntegData.fJacA0  * this->getIntegPointWeightFunction(index);

    VecDouble forcingF(1);
    VecDouble x_ = this->getIntegPointCoordinatesValue(index);
    if (force) force(x_,forcingF);

    for (int i = tshape::NElNodes; i--; ){
        double shapeFi = this->fIntegData.fPhi[i];

        //Viscosity
        double K = 0.;
        for (int l=dim; l--; ) K += this->fIntegData.fDPhiX0(i,l) * du_dx(0,l);

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

template <class tshape>
void ElPoisson<tshape>::ComputeError(VecDouble &errors){
    int index = 0;
    errors.resize(3);
    errors.setZero();
    int DIM = this->Mesh()->Dimension();
    int DEG = this->Mesh()->GetDefaultOrder();

    IntegQuadrature nQuad(DIM,DEG);

    auto exactSol = this->Mesh()->getProblemParameters().getExactSolution();
    if (!exactSol) PanicButton();

    for(int it = 0; it < nQuad.getNumberOfIntegrationPoints(); it++){

        //Defines the integration points adimentional coordinates
        for (int i = DIM; i--; ) this->fIntegData.fAdimCoord[i] = nQuad.PointList(index,i);

        //Returns the quadrature integration weight
        this->fIntegData.fWeight = nQuad.WeightList(index);

        //Computes the jacobian matrix
        this->ComputeJacobian(index);
                    
        this->ComputeSpatialDerivatives();
        
        VecDouble uMEF_(1);
        this->interpolateSolution(index, uMEF_);
        MatrixDouble du_dxMEF(1,DIM);
        this->interpolateSolDerivatives(du_dxMEF);
        
        VecDouble u_(1);
        MatrixDouble gradU(DIM,1);

        VecDouble xna_ = this->getIntegPointCoordinatesValue(index);
        
        exactSol(xna_,u_,gradU);


        //Consider Arlequin weight function
        u_ *= this->getIntegPointWeightFunction(index);
        gradU *= this->getIntegPointWeightFunction(index);
        uMEF_ *= this->getIntegPointWeightFunction(index);
        du_dxMEF *= this->getIntegPointWeightFunction(index);

        //L2 state variable
        errors[0] += (u_[0]-uMEF_[0])*(u_[0]-uMEF_[0]) * this->fIntegData.fWeight * this->fIntegData.fJacA0 ;
        
        //Semi H1 state variable
        for (int m = DIM; m--; ){
            errors[1] += (gradU(m,0)-du_dxMEF(0,m))* (gradU(m,0)-du_dxMEF(0,m)) * this->fIntegData.fWeight * this->fIntegData.fJacA0;
        }

        index++;        
    }; 

    //H1 state variable
    errors[2] = errors[0]+errors[1];
}

template <class tshape>
void ElPoisson<tshape>::ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs){

    for (int i = tshape::NElNodes; i--; ){
        if ((this->Mesh()->NodeVec()[this->getConnectivity()[i]] -> getConstrains(0) == 1) ||
            (this->Mesh()->NodeVec()[this->getConnectivity()[i]] -> getConstrains(0) == 3))  {
            for (int j = tshape::NElNodes; j--; ){
                Stiffness(i,j) = 0.;
                Stiffness(j,i) = 0.;
            };
            Stiffness(i,i) = 1.;
            Rhs[i] = 0.;
        }
    }
    // std::cout<<"Rhs -" << Rhs<<std::endl;

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

template class ElPoisson<ShapePoint>;
template class ElPoisson<ShapeOneD>;
template class ElPoisson<ShapeTriangleLin>;
template class ElPoisson<ShapeTriangleQua>;
template class ElPoisson<ShapeTriangleCub>;
template class ElPoisson<ShapeQuadrilateralLin>;
template class ElPoisson<ShapeTetrahedronLin>;
template class ElPoisson<ShapeTetrahedronQua>;
template class ElPoisson<ShapeTetrahedronCub>;
template class ElPoisson<ShapeHexahedron>;
