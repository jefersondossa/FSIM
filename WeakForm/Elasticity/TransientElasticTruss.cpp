#include "TransientElasticTruss.h"

TransientElasticTruss::TransientElasticTruss(int matid, int dim, double young, double area, double damp, double dens, double dt, TimeIntegScheme tscheme) : ElasticTruss(matid,dim,young,area) {
    fDamping = damp;
    fDensity = dens;
    fTimeStep = dt;
    fIntegScheme = tscheme;
};


void TransientElasticTruss::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    
    MatrixDouble K(Stiffness.rows(),Stiffness.cols());
    K.setZero();
    ElasticTruss::ComputeStiffness(index,data,K);
    
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    int nphi = data.fPhi.size();

    //Compute the mass matrix
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
        Stiffness += K+(1./(fBeta*fTimeStep*fTimeStep) + fGamma*fDamping/(fBeta*fTimeStep)) * Mass;    
        break;
    
    default:
        PanicButton();
        break;
    }
}

void TransientElasticTruss::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){
    VecDouble Residual(Rhs.size());
    Residual.setZero();
    ElasticTruss::ComputeResidual(index,data,Residual);
    
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    int nphi = data.fPhi.size();
    auto vel = data.fDSolDt;
    auto acel = data.fDSolDDt;
    auto disp = data.fSol;

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
    int var = ElasticTruss::VariableIndex(name);

    if (var != -1){
        return var;
    } else {
        if(!strcmp("Velocity",name.c_str()))           return 10;
        if(!strcmp("Acceleration",name.c_str()))       return 11;
    } 

    return -1;
    
};

int TransientElasticTruss::NSolutionVariables(int var) const{
    int nsol = ElasticTruss::NSolutionVariables(var);
    if (nsol != -1){
        return nsol;
    } else {
        if (var == 10 || var == 11) return 3;
    }
    return -1;

};

void TransientElasticTruss::Solution(IntPointData &data, int var, VecDouble &Sol) {

    ElasticTruss::Solution(data,var,Sol);

    //Velocity
    if (var == 10){
        Sol[0] = data.fDSolDt[0];
        Sol[1] = data.fDSolDt[1];
        Sol[2] = 0.;
        return;
    };
    //Acceleration
    if (var == 11){
        Sol[0] = data.fDSolDDt[0];
        Sol[1] = data.fDSolDDt[1];
        Sol[2] = 0.;
        return;
    };
}; 

void TransientElasticTruss::UpdateTimeDerivatives(CompMesh *cmesh){

    switch (fIntegScheme){
    case ENewmark:
        {
            for (int64_t inode = 0; inode < cmesh->NNodes(); inode++){
                //Update Acceleration
                auto acelPrev = cmesh->ConnectVec()[inode]->SolutionDDTime();
                auto velPrev = cmesh->ConnectVec()[inode]->SolutionDTime();
                auto posiPrev = cmesh->ConnectVec()[inode]->PrevSolution();
                auto posi = cmesh->ConnectVec()[inode]->Solution();
                VecDouble acelUpdated(2), velUpdated(2);
                
                //Update Acceleration
                acelUpdated = (posi-posiPrev)/(fBeta*fTimeStep*fTimeStep) - velPrev/(fBeta*fTimeStep) - acelPrev*(1./(2.*fBeta)-1.);

                // //Update Velocity                
                velUpdated = velPrev + fTimeStep*((1.-fGamma)*acelPrev + fGamma*acelUpdated);
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
};