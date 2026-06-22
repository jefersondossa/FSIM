#include "TransientPoisson.h"

//Class constructor
TransientPoisson::TransientPoisson(int matid, int dim, double damp, double dens, double dt, TimeIntegScheme integscheme, int nState) : Poisson(matid,dim,nState) {
    fDamping = damp;
    fDensity = dens;
    fIntegScheme = integscheme;
    fTimeStep = dt;
    this->fType = WeakFormType::kTransientPoisson;
};

void TransientPoisson::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    
    Poisson::ComputeStiffness(index,data,Stiffness);

    //Mass matrix
    int nphi = data.fPhi.size();
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    MatrixDouble Mass(nphi,nphi);
    Mass.setZero();
    for (size_t i = 0; i < nphi; i++){
        for (size_t j = 0; j < nphi; j++){
            Mass(i,j) += data.fPhi[i] * data.fPhi[j] *  WJ * fDensity;
        }
    }

    switch (fIntegScheme)
    {
    case ENewmark:
        Stiffness += (1./(fBeta*fTimeStep*fTimeStep)) * Mass;    
        // Stiffness += (1./(fBeta*fTimeStep*fTimeStep) + fGamma*fDamping/(fBeta*fTimeStep)) * Mass;    
        break;
    
    default:
        PanicButton();
        break;
    }
}

void TransientPoisson::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    Poisson::ComputeResidual(index,data,Rhs);

    auto vel = data.fDSolDt;
    auto acel = data.fDSolDDt;
    auto disp = data.fSol;
    auto dispPrev = data.fSolPrev;
    int nphi = data.fPhi.size();
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];

    for (size_t i = 0; i < nphi; i++){
        Rhs[i] += (disp[0]/(fBeta * fTimeStep * fTimeStep) + 
                   vel[0]/(fBeta*fTimeStep) + 
                   (1./(2.*fBeta)-1.) * acel[0]) * data.fPhi[i] * fDensity * WJ;
    }

};

void TransientPoisson::ComputeError(IntPointData &data, VecDouble &errors){

    PanicButton();
}


int TransientPoisson::VariableIndex(const std::string &name) const{
    
    return Poisson::VariableIndex(name);
    
};

int TransientPoisson::NSolutionVariables(int var) const{
    return Poisson::NSolutionVariables(var);
};

void TransientPoisson::Solution(IntPointData &data, int var, VecDouble &Sol) {
    Poisson::Solution(data,var,Sol);

}; 

void TransientPoisson::UpdateTimeDerivatives(CompMesh *cmesh){
    switch (fIntegScheme){
    case ENewmark:
        {
            for (int inode = 0; inode < cmesh->NConnects(); inode++){
                //Update Acceleration
                auto acelPrev = cmesh->ConnectVec()[inode]->SolutionDDTime();
                auto velPrev = cmesh->ConnectVec()[inode]->SolutionDTime();
                auto dispPrev = cmesh->ConnectVec()[inode]->PrevSolution();
                auto disp = cmesh->ConnectVec()[inode]->Solution();
                VecDouble acelUpdated(2), velUpdated(2);
                acelUpdated = (disp-dispPrev)/(fBeta*fTimeStep*fTimeStep) -
                                velPrev/(fBeta*fTimeStep) -
                                (1./(2.*fBeta) - 1.) * acelPrev; 

                //Update Velocity
                velUpdated = velPrev + (1.-fGamma)*fTimeStep*acelPrev + fGamma*fTimeStep*acelUpdated;

                cmesh->ConnectVec()[inode]->SetDSolutionDTime(0,velUpdated[0]);
                cmesh->ConnectVec()[inode]->SetDSolutionDDTime(0,acelUpdated[0]);
            }           
        }
        break;
    case EGeneralizedAlpha:
        PanicButton();
        // Nesse caso vai ser complicado porque o domínio /gradiente da função mudança
        // de configuração deve ser interpolada entre os passos anterior e atual
        /* code */
        break;
    
    default:
        PanicButton();
        break;
    }
}