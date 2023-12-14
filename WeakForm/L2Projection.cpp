#include "L2Projection.h"

L2Projection::L2Projection(int matid, int dim, int bctype, MatrixDouble &val1, VecDouble &val2) : WeakForm() {
    this->fDimension = dim;
    this->fMatId = matid;
    fNState = val2.size();
    BCType = bctype;
    BCVal1 = val1;
    BCVal2 = val2;
};

void L2Projection::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    if (!data.fNeedsSol){
        data.fNeedsSol = true;
        data.fSol.resize(fNState);
    }
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index] * WeakForm::fBigNumber;
    double nphi = data.fPhi.size();

    switch (BCType)
    {
    case 0: // Dirichlet
        for (int i = 0; i < nphi; i++){
            for (int j = 0; j < nphi; j++){
                for (int istate = 0; istate < fNState; istate++){
                    Stiffness(fNState*i+istate,fNState*j+istate) +=  WJ * data.fPhi[i] * data.fPhi[j];
                }
            }
        }
        break;
    case 1:
        break;
    case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
        for(int i = 0 ; i < nphi; i++) {
//                ef(nstate*in+0,0) += BIGNUMBER * (0. - data.sol[0][0]) * v2[0] * phi(in,0) * weight;
//                ef(nstate*in+1,0) += BIGNUMBER * (0. - data.sol[0][1]) * v2[1] * phi(in,0) * weight;
            for (int j = 0 ; j < nphi; j++) {
                for (int istate = 0; istate < fNState; istate++){
                    Stiffness(fNState*i+istate,fNState*j+istate) += WJ * data.fPhi[i] * data.fPhi[j] * BCVal2[istate];
                }
            }//jn
        }//in
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
    case 0: // Dirichlet in all state variables
        for (int i = 0; i < nphi; i++){
            for (int istate = 0; istate < fNState; istate++){
                Rhs(fNState*i+istate) +=  WeakForm::fBigNumber * WJ * data.fPhi[i] * result[istate];
            }
        }
        break;
    case 1: // Neumann in all state variables
    {
        for (int i = 0; i < nphi; i++){
            for (int istate = 0; istate < fNState; istate++){
                Rhs(fNState*i+istate) +=  WJ * data.fPhi[i] * result[istate];
            }
        }
    }
        break;
    case 3:
        break;
    default:
        std::cout << "BC Type not implemented \n" ;
        PanicButton();
        break;
    }
};

void L2Projection::ComputeError(IntPointData &data, VecDouble &errors){

}



