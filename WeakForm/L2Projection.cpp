#include "L2Projection.h"

L2Projection::L2Projection(int matid, int dim, BoundaryConditionType bctype, MatrixDouble &val1, VecDouble &val2) : WeakForm() {
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
    
    MatrixDouble deriv(data.fX.size(), fNState);
    deriv.setZero();
    VecDouble result = BCVal2;
    
    if(this->fExactSol)
    {
        fExactSol(data.fX, result, deriv);
    }
    // result -= data.fSol;

    switch (BCType)
    {
    case BoundaryConditionType::kDirichlet: // Dirichlet
        for (int i = 0; i < nphi; i++){
            for (int j = 0; j < nphi; j++){
                for (int istate = 0; istate < fNState; istate++){
                    Stiffness(fNState*i+istate,fNState*j+istate) +=  WJ * data.fPhi[i] * data.fPhi[j];
                }
            }
        }
        break;
    case BoundaryConditionType::kNeumann:
        break;
    case BoundaryConditionType::kDirectionalHomogeneousDirichlet: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
        for(int i = 0 ; i < nphi; i++) {
            for (int j = 0 ; j < nphi; j++) {
                for (int istate = 0; istate < fNState; istate++){
                    Stiffness(fNState*i+istate,fNState*i+istate) += WJ * data.fPhi[i] * data.fPhi[j] * BCVal2[istate];
                }
            }//jn
        }//in
        break;
    case BoundaryConditionType::kDirectionalNonHomogeneousDirichlet: // Directional NonHomogeneous Dirichlet - displacement is set to null in the non-null vector component direction
        for(int i = 0 ; i < nphi; i++) {
            for (int j = 0 ; j < nphi; j++) {
                for (int istate = 0; istate < fNState; istate++){
                    if (fabs(BCVal2[istate]) > 0 && fabs(result[istate])>1.e-10)
                    Stiffness(fNState*i+istate,fNState*i+istate) += WJ * data.fPhi[i] * data.fPhi[j];
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
    
    auto force = fForceFunction;
    VecDouble forcingF(3);
    forcingF.setZero();
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    auto forceT = fForceFunctionTransient;
    VecDouble forcingFT(3);
    forcingFT.setZero();
    double time = fTimeInstant * fTimeStep;
    if (forceT) forceT(x_, time, forcingFT);

    for (int i = 0; i < fNState; i++){
        result[i] += forcingF[i] + forcingFT[i];
        if (this->BCType == BoundaryConditionType::kNeumann){
            BCVal2[i] -= forcingF[i] + forcingFT[i]; // Add time derivative contribution
        }
    }
    

    switch (BCType)
    {
    case BoundaryConditionType::kDirichlet: // Dirichlet in all state variables
        for (int i = 0; i < nphi; i++){
            for (int istate = 0; istate < fNState; istate++){
                Rhs(fNState*i+istate) +=  WeakForm::fBigNumber * WJ * data.fPhi[i] * result[istate];
            }
        }
        break;
    case BoundaryConditionType::kNeumann: // Neumann in all state variables
    {
        for (int i = 0; i < nphi; i++){
            for (int istate = 0; istate < fNState; istate++){
                Rhs(fNState*i+istate) +=  WJ * data.fPhi[i] * (BCVal2[istate]);
            }
        }
    }
        break;
    case BoundaryConditionType::kDirectionalHomogeneousDirichlet:
        break;
    
    case BoundaryConditionType::kDirectionalNonHomogeneousDirichlet:
        for (int i = 0; i < nphi; i++){
            for (int istate = 0; istate < fNState; istate++){
                if (fabs(BCVal2[istate]) > 0 && fabs(result[istate])>1.e-10){
                    Rhs(fNState*i+istate) +=  WeakForm::fBigNumber * WJ * data.fPhi[i] * (result[istate]);
                }
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

}



int L2Projection::VariableIndex(const std::string &name) const{
    
    if(!strcmp("Solution",name.c_str()))        return 1;
    if(!strcmp("DerivativeX",name.c_str()))      return 2;
    if(!strcmp("DerivativeY",name.c_str()))      return 3;
    if(!strcmp("DerivativeZ",name.c_str()))      return 4;
    if(!strcmp("Material",name.c_str()))         return 100;
    
    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int L2Projection::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
    case 2:
    case 3:
    case 4:
    case 100:
        return 1;

    default:
        PanicButton();
        return -1;
    }
};

void L2Projection::Solution(IntPointData &data, int var, VecDouble &Sol) {

    if (var == 100){
        Sol[0] = fMatId;
        return;
    };

    //Solution
    if (var == 1){
        for (int i = 0; i < fNState; i++){
            Sol[i] = data.fSol[i];
        }
        return;
    };

    //Derivative X
    if (var == 2){
        for (int i = 0; i < fNState; i++){
            Sol[i] = data.fDSolDx(0,i);
        }
        return;
    };
    //Derivative Y
    if (var == 3){
        for (int i = 0; i < fNState; i++){
            Sol[i] = data.fDSolDx(1,i);
        }
        return;
    };
    //Derivative Z
    if (var == 4){
        for (int i = 0; i < fNState; i++){
            Sol[i] = data.fDSolDx(2,i);
        }
        return;
    };

}; 