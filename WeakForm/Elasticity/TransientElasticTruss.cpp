#include "TransientElasticTruss.h"

TransientElasticTruss::TransientElasticTruss(int matid, int dim, double young, double area, double damp, double dens, double dt, TimeIntegScheme tscheme) : ElasticTruss(matid,dim,young,area) {
    fDamping = damp;
    fDensity = dens;
    fTimeStep = dt;
};


void TransientElasticTruss::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    ElasticTruss::ComputeStiffness(index,data,Stiffness);
    
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    int nphi = data.fPhi.size();
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
        PanicButton();
        break;
    }
}

void TransientElasticTruss::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){
    ElasticTruss::ComputeResidual(index,data,Rhs);
    
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    int nphi = data.fPhi.size();
    auto vel = data.fDSolDt;
    auto acel = data.fDSolDDt;
    auto disp = data.fSol;
    auto dispPrev = data.fSolPrev;

    for (size_t i = 0; i < nphi; i++){
        Rhs[2*i  ] += (disp[0]/(fBeta * fTimeStep * fTimeStep) + 
                       vel[0]/(fBeta*fTimeStep) + 
                       (1./(2.*fBeta)-1.) * acel[0]) * data.fPhi[i] * fDensity * WJ +
                       (disp[0]*fGamma/(fBeta * fTimeStep) +
                       vel[0] * (fGamma/fBeta - 1.) +
                       acel[0] * fTimeStep * (fGamma/(2.*fBeta)-1.)) * data.fPhi[i] * fDensity * fDamping * WJ;
        Rhs[2*i+1] += (disp[1]/(fBeta * fTimeStep * fTimeStep) + 
                       vel[1]/(fBeta*fTimeStep) + 
                       (1./(2.*fBeta)-1.) * acel[1]) * data.fPhi[i] * fDensity * WJ +
                       (disp[1]*fGamma/(fBeta * fTimeStep) +
                       vel[1] * (fGamma/fBeta - 1.) +
                       acel[1] * fTimeStep * (fGamma/(2.*fBeta)-1.)) * data.fPhi[i] * fDensity * fDamping * WJ;
    }
};

void TransientElasticTruss::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();

}

int TransientElasticTruss::VariableIndex(const std::string &name) const{
    return ElasticTruss::VariableIndex(name);
    
};

int TransientElasticTruss::NSolutionVariables(int var) const{
    return ElasticTruss::NSolutionVariables(var);

};

void TransientElasticTruss::Solution(IntPointData &data, int var, VecDouble &Sol) {

    ElasticTruss::Solution(data,var,Sol);

}; 

void TransientElasticTruss::UpdateTimeDerivatives(CompMesh *cmesh){

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