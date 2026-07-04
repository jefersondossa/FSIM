#include "TransientPositionalFrame2D.h"

TransientPositionalFrame2D::TransientPositionalFrame2D(int matid, REAL young, REAL base, REAL height, REAL damp, REAL dens, REAL dt, TimeIntegScheme tscheme) : PositionalFrame2D(matid,young,base,height) {    
    fDamping = damp;
    fDensity = dens;
    fTimeStep = dt;
    this->fType = WeakFormType::kTransientPositionalFrame2D;
};

void TransientPositionalFrame2D::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    
    PositionalFrame2D::ComputeStiffness(index,data,Stiffness);

    int nphi = data.fPhi.size();
    REAL WJ = data.fWeight * data.fJacA0;
    
    REAL initLenght = 2.*data.fJacA0;
    MatrixDouble Mass(3*nphi,3*nphi);
    Mass.setZero();
    for (size_t i = 0; i < nphi; i++){
        for (size_t j = 0; j < nphi; j++){
            Mass(3*i  ,3*j  ) += data.fPhi[i] * data.fPhi[j] *  WJ * fDensity * fDepth;
            Mass(3*i+1,3*j+1) += data.fPhi[i] * data.fPhi[j] *  WJ * fDensity * fDepth;
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
    REAL WJ = data.fWeight * data.fJacA0 * fDepth * fDensity;
    
    auto vel = data.fDSolDt;
    auto acel = data.fDSolDDt;
    auto posi = data.fSol;
    auto posiPrev = data.fSolPrev;
    VecDouble Rhs2(3*nphi);
    Rhs2.setZero();

    auto force = fForceFunctionTransient;
    VecDouble forcingF(3);
    forcingF.setZero();
    VecDouble x_ = data.fX;
    REAL time = fTimeInstant * fTimeStep;
    if (force) force(x_, time, forcingF);

    REAL a0 = 1./(fBeta * fTimeStep * fTimeStep);
    REAL a1 = fGamma / (fBeta * fTimeStep);
    REAL a2 = 1./(fTimeStep * fBeta);
    REAL a3 = 1. / (2. * fBeta);
    REAL a4 = fGamma / fBeta;
    REAL a5 = fTimeStep * (fGamma / (2. * fBeta) - 1);
    
    for (size_t i = 0; i < nphi; i++){
        Rhs[3*i  ] += forcingF[0] * data.fPhi[i] * WJ
                      -((posi[0]-posiPrev[0]) * (1/(fBeta * fTimeStep * fTimeStep) + fDamping * fGamma / (fBeta * fTimeStep))) * data.fPhi[i] * WJ
                      + data.fPhi[i] * WJ * (a2 * vel[0] + a3 * acel[0] ) +  data.fPhi[i] * WJ*fDamping * (a4 * vel[0] + a5 * acel[0]);
                    
        Rhs[3*i+1] += forcingF[1] * data.fPhi[i] * WJ
                      -((posi[1]-posiPrev[1]) * (1/(fBeta * fTimeStep * fTimeStep) + fDamping * fGamma / (fBeta * fTimeStep))) * data.fPhi[i] * WJ
                      + data.fPhi[i] * WJ * (a2 * vel[1] + a3 * acel[1] ) +  data.fPhi[i] * WJ*fDamping * (a4 * vel[1] + a5 * acel[1]);
    }
    // std::cout << "rhs " << Rhs[3] << ", Rhs2 " << Rhs2[3] << std::endl; 
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
            for (int inode = 0; inode < cmesh->NConnects(); inode++){
                //Update Acceleration
                auto acelPrev = cmesh->ConnectVec()[inode]->SolutionDDTime();
                auto velPrev = cmesh->ConnectVec()[inode]->SolutionDTime();
                auto posiPrev = cmesh->ConnectVec()[inode]->PrevSolution();
                auto posi = cmesh->ConnectVec()[inode]->Solution();
               
               
                VecDouble acelUpdated(2), velUpdated(2);

                // //Update Velocity    
                velUpdated = velPrev + (fGamma/(fBeta * fTimeStep)) * posi
                             - (fGamma / fBeta) * velPrev 
                             - fTimeStep * (fGamma / (2 * fBeta) - 1) * acelPrev;
                
                //Update Acceleration
                acelUpdated = acelPrev + (1 / (fBeta * fTimeStep * fTimeStep)) * posi
                              - 1 / (fBeta * fTimeStep) * velPrev
                              - 1 / (2 * fBeta) * acelPrev;

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