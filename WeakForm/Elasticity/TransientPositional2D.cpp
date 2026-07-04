#include "TransientPositional2D.h"

TransientPositional2D::TransientPositional2D(int matid, REAL young, REAL poisson, bool planes, REAL damp, REAL dens, REAL dt, TimeIntegScheme integscheme) : 
                        ElasticityPositional2D(matid,young,poisson,planes) {
    fDamping = damp;
    fDensity = dens;
    fTimeStep = dt;
    fIntegScheme = integscheme;
    this->fType = WeakFormType::kTransientPositional2D;
};

void TransientPositional2D::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    ElasticityPositional2D::ComputeStiffness(index,data,Stiffness);

    int nphi = data.fPhi.size();
    REAL WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    
    REAL initLenght = 2.*data.fJacA0;
    MatrixDouble Mass(2*nphi,2*nphi);
    Mass.setZero();
    for (size_t i = 0; i < nphi; i++){
        for (size_t j = 0; j < nphi; j++){
            Mass(2*i  ,2*j  ) += data.fPhi[i] * data.fPhi[j] *  WJ * fDensity;
            Mass(2*i+1,2*j+1) += data.fPhi[i] * data.fPhi[j] *  WJ * fDensity;
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

void TransientPositional2D::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    ElasticityPositional2D::ComputeResidual(index,data,Rhs);

    int nphi = data.fPhi.size();
    REAL WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    
    auto vel = data.fDSolDt;
    auto acel = data.fDSolDDt;
    auto disp = data.fSol;
    auto dispPrev = data.fSolPrev;

    auto qs = dispPrev/(fBeta*fTimeStep*fTimeStep) + vel/(fBeta*fTimeStep) +
                    (1./(2.*fBeta) - 1.) * acel;
                    
    for (size_t i = 0; i < nphi; i++){
        Rhs[2*i  ] -= (disp[0]/(fBeta * fTimeStep * fTimeStep) -qs[0]) * data.fPhi[i] * fDensity * WJ;
        Rhs[2*i+1] -= (disp[1]/(fBeta * fTimeStep * fTimeStep) -qs[1]) * data.fPhi[i] * fDensity * WJ;
    }
};

void TransientPositional2D::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();

}

int TransientPositional2D::VariableIndex(const std::string &name) const{
    
    return ElasticityPositional2D::VariableIndex(name);
};

int TransientPositional2D::NSolutionVariables(int var) const{
    return ElasticityPositional2D::NSolutionVariables(var);
};

void TransientPositional2D::Solution(IntPointData &data, int var, VecDouble &Sol) {

    ElasticityPositional2D::Solution(data,var,Sol);

}; 


void TransientPositional2D::UpdateTimeDerivatives(CompMesh *cmesh){

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
                cmesh->ConnectVec()[inode]->SetDSolutionDTime(1,velUpdated[1]);

                cmesh->ConnectVec()[inode]->SetDSolutionDDTime(0,acelUpdated[0]);
                cmesh->ConnectVec()[inode]->SetDSolutionDDTime(1,acelUpdated[1]);
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

}

