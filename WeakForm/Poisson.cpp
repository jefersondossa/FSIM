#include "Poisson.h"

//Class constructor
Poisson::Poisson(int matid, int dim, int nState) : WeakForm() {
    fMatId = matid;
    fDimension = dim;
    fNState = nState;
};

void Poisson::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    if (!data.fNeedsDSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
    }

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    int nphi = data.fPhi.size();

    for (int i = nphi; i-- ; ){
        for (int j = nphi; j-- ; ){            
            for (int k = fDimension; k--;  ){
                Stiffness(i,j) += data.fDPhiX0(i,k) * data.fDPhiX0(j,k) * WJ;
            }
        };
    };
}

void Poisson::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    auto force = fForceFunction;
    int nphi = data.fPhi.size();

    double WJ = data.fWeight * data.fJacA0  * data.fWeightFunction[index];

    VecDouble forcingF(1);
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    for (int i = nphi; i--; ){
        double shapeFi = data.fPhi[i];

        //Matrix residual
        double K = 0.;
        for (int l=fDimension; l--; ) K += data.fDPhiX0(i,l) * data.fDSolDx(0,l);

        //Source term
        double F = (forcingF[0]) * shapeFi;

        Rhs[i] += (-K + F) * WJ;
    };
};

void Poisson::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Need refactor \n" << std::endl;
    PanicButton();
    // int index = 0;
    // errors.resize(3);
    // errors.setZero();
    // int DIM = this->Mesh()->Dimension();
    // int DEG = this->Mesh()->GetDefaultOrder();

    // IntegQuadrature nQuad(DIM,DEG);

    // auto exactSol = this->Mesh()->getProblemParameters().getExactSolution();
    // if (!exactSol) PanicButton();

    // for(int it = 0; it < nQuad.getNumberOfIntegrationPoints(); it++){

    //     //Defines the integration points adimentional coordinates
    //     for (int i = DIM; i--; ) data.fAdimCoord[i] = nQuad.PointList(index,i);

    //     //Returns the quadrature integration weight
    //     data.fWeight = nQuad.WeightList(index);

    //     //Computes the jacobian matrix
    //     this->ComputeJacobian(index);
                    
    //     this->ComputeSpatialDerivatives();
        
    //     VecDouble uMEF_(1);
    //     this->interpolateSolution(index, uMEF_);
    //     MatrixDouble du_dxMEF(1,DIM);
    //     this->interpolateSolDerivatives(du_dxMEF);
        
    //     VecDouble u_(1);
    //     MatrixDouble gradU(DIM,1);

    //     VecDouble xna_ = this->getIntegPointCoordinatesValue(index);
        
    //     exactSol(xna_,u_,gradU);


    //     //Consider Arlequin weight function
    //     u_ *= data.fWeightFunction[index];
    //     gradU *= data.fWeightFunction[index];
    //     uMEF_ *= data.fWeightFunction[index];
    //     du_dxMEF *= data.fWeightFunction[index];

    //     //L2 state variable
    //     errors[0] += (u_[0]-uMEF_[0])*(u_[0]-uMEF_[0]) * data.fWeight * data.fJacA0 ;
        
    //     //Semi H1 state variable
    //     for (int m = DIM; m--; ){
    //         errors[1] += (gradU(m,0)-du_dxMEF(0,m))* (gradU(m,0)-du_dxMEF(0,m)) * data.fWeight * data.fJacA0;
    //     }

    //     index++;        
    // }; 

    // //H1 state variable
    // errors[2] = errors[0]+errors[1];
}


