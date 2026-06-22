#include "TransientPositionalTruss.h"

TransientPositionalTruss::TransientPositionalTruss(int matid, int dim, double young, double area, double damp, double dens, double dt, TimeIntegScheme tscheme) : PositionalTruss(matid,dim,young,area) {    
    fDamping = damp;
    fDensity = dens;
    fTimeStep = dt;
    this->fType = WeakFormType::kTransientPositionalTruss;
};

void TransientPositionalTruss::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    
    PositionalTruss::ComputeStiffness(index,data,Stiffness);

    int nphi = data.fPhi.size();
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    
    double initLenght = 2.*data.fJacA0;
    MatrixDouble Mass(2*nphi,2*nphi);
    Mass.setZero();
    for (size_t i = 0; i < nphi; i++){
        for (size_t j = 0; j < nphi; j++){
            Mass(2*i  ,2*j  ) += data.fPhi[i] * data.fPhi[j] *  WJ * fDensity * fArea;
            Mass(2*i+1,2*j+1) += data.fPhi[i] * data.fPhi[j] *  WJ * fDensity * fArea;
        }
    }

    switch (fIntegScheme)
    {
    case ENewmark:
        Stiffness += (1./(fBeta*fTimeStep*fTimeStep) + fGamma*fDamping/(fBeta*fTimeStep)) * Mass;    
        break;
    
    default:
        break;
    }

}

void TransientPositionalTruss::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    PositionalTruss::ComputeResidual(index,data,Rhs);
    int nphi = data.fPhi.size();
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index] * fArea * fDensity;
    
    auto vel = data.fDSolDt;
    auto acel = data.fDSolDDt;
    auto posi = data.fSol;
    auto posiPrev = data.fSolPrev;

    VecDouble qs = posiPrev/(fBeta*fTimeStep*fTimeStep) + vel/(fBeta*fTimeStep) +
                    (1./(2.*fBeta) - 1.) * acel;
    for (size_t i = 0; i < nphi; i++){
        Rhs[2*i  ] -= (posi[0]/(fBeta * fTimeStep * fTimeStep) - qs[0]) * data.fPhi[i] * WJ;
        Rhs[2*i+1] -= (posi[1]/(fBeta * fTimeStep * fTimeStep) - qs[1]) * data.fPhi[i] * WJ;
    }
    
};

void TransientPositionalTruss::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();

}


int TransientPositionalTruss::VariableIndex(const std::string &name) const{
    
    return PositionalTruss::VariableIndex(name);
};

int TransientPositionalTruss::NSolutionVariables(int var) const{
    return PositionalTruss::NSolutionVariables(var);
};

void TransientPositionalTruss::Solution(IntPointData &data, int var, VecDouble &Sol) {

    PositionalTruss::Solution(data,var,Sol);

}; 


void TransientPositionalTruss::UpdateTimeDerivatives(CompMesh *cmesh){

    switch (fIntegScheme){
    case ENewmark:
        {
            for (int inode = 0; inode < cmesh->NConnects(); inode++){
                //Update Acceleration
                auto acelPrev = cmesh->ConnectVec()[inode]->SolutionDDTime();
                auto velPrev = cmesh->ConnectVec()[inode]->SolutionDTime();
                auto posiPrev = cmesh->ConnectVec()[inode]->PrevSolution();
                auto posi = cmesh->ConnectVec()[inode]->Solution();
                VecDouble acelUpdated(2), velUpdated(2);
                auto qs = posiPrev/(fBeta*fTimeStep*fTimeStep) + velPrev/(fBeta*fTimeStep) +
                                (1./(2.*fBeta) - 1.) * acelPrev;
                auto rs = velPrev + (1.-fGamma)*fTimeStep*acelPrev;
                //Update Acceleration
                acelUpdated = posi/(fBeta*fTimeStep*fTimeStep) - qs;

                // //Update Velocity                
                velUpdated = posi*fGamma/(fBeta*fTimeStep) + rs -fGamma*fTimeStep*qs;

                cmesh->ConnectVec()[inode]->SetDSolutionDTime(0,velUpdated[0]);
                cmesh->ConnectVec()[inode]->SetDSolutionDTime(1,velUpdated[1]);

                cmesh->ConnectVec()[inode]->SetDSolutionDDTime(0,acelUpdated[0]);
                cmesh->ConnectVec()[inode]->SetDSolutionDDTime(1,acelUpdated[1]);
                // std::cout << "Acel = " << acelUpdated[0] << std::endl;
                // std::cout << "Vel = " << velUpdated[0] << std::endl;
            }           
        }
        break;
    case EGeneralizedAlpha:
        // Nesse caso vai ser complicado porque o domínio /gradiente da função mudança
        // de configuração deve ser interpolada entre os passos anterior e atual
        /* code */
        break;
    
    default:
        PanicButton();
        break;
    }
};