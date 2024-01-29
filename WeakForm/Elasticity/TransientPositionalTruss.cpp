#include "TransientPositionalTruss.h"

TransientPositionalTruss::TransientPositionalTruss(int matid, int dim, double young, double area, double damp, double dens, double dt, TimeIntegScheme tscheme) : PositionalTruss(matid,dim,young,area) {
    std::cout << "This material is not working properly. Please debug it" << std::endl;
    
    fDamping = damp;
    fDensity = dens;
    fTimeStep = dt;
};

void TransientPositionalTruss::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    
    PositionalTruss::ComputeStiffness(index,data,Stiffness);

    int nphi = data.fPhi.size();
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    
    double initLenght = 2.*data.fJacA0;
    MatrixDouble Mass(2*nphi,2*nphi);
    Mass.setIdentity();

    Mass *= initLenght*fArea*fDensity*0.5;

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
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    
    auto vel = data.fDSolDt;
    auto acel = data.fDSolDDt;
    auto disp = data.fSol;
    auto dispPrev = data.fSolPrev;

    auto qs = dispPrev/(fBeta*fTimeStep*fTimeStep) + vel/(fBeta*fTimeStep) +
                    (1./(2.*fBeta) - 1.) * acel;
    auto rs = vel + (1.-fGamma)*fTimeStep*acel;

    for (size_t i = 0; i < nphi; i++){
        Rhs[2*i  ] += (disp[0]/(fBeta * fTimeStep * fTimeStep) -qs[0]) * data.fPhi[i] * fDensity * WJ;
        Rhs[2*i+1] += (disp[1]/(fBeta * fTimeStep * fTimeStep) -qs[1]) * data.fPhi[i] * fDensity * WJ;
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
            for (int64_t inode = 0; inode < cmesh->NNodes(); inode++){
                //Update Acceleration
                auto acelPrev = cmesh->NodeVec()[inode]->SolutionDDTime();
                auto velPrev = cmesh->NodeVec()[inode]->SolutionDTime();
                auto dispPrev = cmesh->NodeVec()[inode]->PrevSolution();
                auto disp = cmesh->NodeVec()[inode]->Solution();
                VecDouble acelUpdated(2), velUpdated(2);
                auto qs = dispPrev/(fBeta*fTimeStep*fTimeStep) + velPrev/(fBeta*fTimeStep) +
                                (1./(2.*fBeta) - 1.) * acelPrev;
                auto rs = velPrev + (1.-fGamma)*fTimeStep*acelPrev;

                // acelUpdated = (disp-dispPrev)/(fBeta*fTimeStep*fTimeStep) -
                //                 velPrev/(fBeta*fTimeStep) -
                //                 (1./(2.*fBeta) - 1.) * acelPrev; 
                acelUpdated = disp/(fBeta*fTimeStep*fTimeStep) - qs;

                //Update Velocity
                // velUpdated = velPrev + (1.-fGamma)*fTimeStep*acelPrev + fGamma*fTimeStep*acelUpdated;
                
                velUpdated = disp*fGamma/(fBeta*fTimeStep) + rs -fGamma*fTimeStep*qs;

                cmesh->NodeVec()[inode]->SetDSolutionDTime(0,velUpdated[0]);
                cmesh->NodeVec()[inode]->SetDSolutionDTime(1,velUpdated[1]);
                
                cmesh->NodeVec()[inode]->SetDSolutionDDTime(0,acelUpdated[0]);
                cmesh->NodeVec()[inode]->SetDSolutionDDTime(1,acelUpdated[1]);
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