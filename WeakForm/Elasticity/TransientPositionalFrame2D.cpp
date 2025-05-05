#include "TransientPositionalFrame2D.h"

TransientPositionalFrame2D::TransientPositionalFrame2D(int matid, double young, double base, double height, double damp, double dens, double dt, TimeIntegScheme tscheme) : PositionalFrame2D(matid,young,base,height) {    
    fDamping = damp;
    fDensity = dens;
    fTimeStep = dt;
};

void TransientPositionalFrame2D::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    
    PositionalFrame2D::ComputeStiffness(index,data,Stiffness);

    int nphi = data.fPhi.size();
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    
    double initLenght = 2.*data.fJacA0;
    MatrixDouble Mass(2*nphi,2*nphi);
    Mass.setZero();
    for (size_t i = 0; i < nphi; i++){
        for (size_t j = 0; j < nphi; j++){
            Mass(2*i  ,2*j  ) += data.fPhi[i] * data.fPhi[j] *  WJ * fDensity * fDepth;
            Mass(2*i+1,2*j+1) += data.fPhi[i] * data.fPhi[j] *  WJ * fDensity * fDepth;
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

void TransientPositionalFrame2D::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    PositionalFrame2D::ComputeResidual(index,data,Rhs);
    int nphi = data.fPhi.size();
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index] * fDepth * fDensity;
    
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

void TransientPositionalFrame2D::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();

}


int TransientPositionalFrame2D::VariableIndex(const std::string &name) const{
    
    return PositionalFrame2D::VariableIndex(name);
};

int TransientPositionalFrame2D::NSolutionVariables(int var) const{
    return PositionalFrame2D::NSolutionVariables(var);
};

void TransientPositionalFrame2D::Solution(IntPointData &data, int var, VecDouble &Sol) {

    PositionalFrame2D::Solution(data,var,Sol);

}; 


void TransientPositionalFrame2D::UpdateTimeDerivatives(CompMesh *cmesh){

    switch (fIntegScheme){
    case ENewmark:
        {
            for (int64_t inode = 0; inode < cmesh->NNodes(); inode++){
                //Update Acceleration
                auto acelPrev = cmesh->NodeVec()[inode]->SolutionDDTime();
                auto velPrev = cmesh->NodeVec()[inode]->SolutionDTime();
                auto posiPrev = cmesh->NodeVec()[inode]->PrevSolution();
                auto posi = cmesh->NodeVec()[inode]->Solution();
                VecDouble acelUpdated(2), velUpdated(2);
                auto qs = posiPrev/(fBeta*fTimeStep*fTimeStep) + velPrev/(fBeta*fTimeStep) +
                                (1./(2.*fBeta) - 1.) * acelPrev;
                auto rs = velPrev + (1.-fGamma)*fTimeStep*acelPrev;
                //Update Acceleration
                acelUpdated = posi/(fBeta*fTimeStep*fTimeStep) - qs;

                // //Update Velocity                
                velUpdated = posi*fGamma/(fBeta*fTimeStep) + rs -fGamma*fTimeStep*qs;

                cmesh->NodeVec()[inode]->SetDSolutionDTime(0,velUpdated[0]);
                cmesh->NodeVec()[inode]->SetDSolutionDTime(1,velUpdated[1]);
                
                cmesh->NodeVec()[inode]->SetDSolutionDDTime(0,acelUpdated[0]);
                cmesh->NodeVec()[inode]->SetDSolutionDDTime(1,acelUpdated[1]);
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