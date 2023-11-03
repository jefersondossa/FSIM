#include "L2Projection.h"

void L2Projection::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    data.fNeedsSol = true;
    data.fSol.resize(fNState);
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index] * WeakForm::fBigNumber;
    double nphi = data.fPhi.size();

    switch (BCType)
    {
    case 0: // Dirichlet
        for (int i = 0; i < nphi; i++){
            for (int j = 0; j < nphi; j++){
                for (int istate = 0; istate < fNState; istate++){
                    for (int jstate = 0; jstate < fNState; jstate++){
                        Stiffness(fNState*i+istate,fNState*j+jstate) +=  WJ * data.fPhi[i] * data.fPhi[j];
                    }
                }
            }
        }
        break;
    case 1:
        break;
    
    default:
        std::cout << "BC Type not implemented \n" ;
        PanicButton();
        break;
    }
   
}

void L2Projection::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    double nphi = data.fPhi.size();

    MatrixDouble deriv(data.fX.size(), fNState);
    deriv.setZero();
    VecDouble result = BCVal2;
    
    if(this->fExactSol)
    {
        fExactSol(data.fX, result, deriv);
    }

    result -= data.fSol;
    switch (BCType)
    {
    case 0: // Dirichlet
        for (int i = 0; i < nphi; i++){
            for (int istate = 0; istate < fNState; istate++){
                Rhs(fNState*i+istate) +=  WeakForm::fBigNumber * WJ * data.fPhi[i] * result[istate];
            }
        }
        break;
    case 1: // Neumann
        for (int i = 0; i < nphi; i++){
            for (int istate = 0; istate < fNState; istate++){
                Rhs(fNState*i+istate) +=  WJ * data.fPhi[i] * result[istate];
            }
        }
        break;
    
    default:
        std::cout << "BC Type not implemented \n" ;
        PanicButton();
        break;
    }
};

void L2Projection::ComputeError(IntPointData &data, VecDouble &errors){
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



