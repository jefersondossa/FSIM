#include "TransientElasticTruss.h"

TransientElasticTruss::TransientElasticTruss(int matid, int dim, double young, double area, double damp, double dens, double dt, TimeIntegScheme tscheme) : ElasticTruss(matid,dim,young,area) {
    std::cout << "This material is not working properly. Please debug it" << std::endl;
    
    fDamping = damp;
    fDensity = dens;
    fTimeStep = dt;
};


void TransientElasticTruss::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    ElasticTruss::ComputeStiffness(index,data,Stiffness);
    
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    int nphi = data.fPhi.size();
    MatrixDouble rotation(fDimension*nphi,fDimension*nphi);
    MatrixDouble matB(fDimension,fDimension*nphi);
    rotation.setZero();
    matB.setZero();
    double cosa = data.fAxes0(0,0) / data.fJacA0;
    double sina = data.fAxes0(1,0) / data.fJacA0;
    double check = sina*sina+cosa*cosa;
    for (int j = 0; j < nphi; j++){
        matB(0,fDimension*j) = data.fPhi[j];
        
        rotation(2*j  ,2*j  ) = cosa;
        rotation(2*j+1,2*j  ) = sina;
        rotation(2*j  ,2*j+1) = -sina;
        rotation(2*j+1,2*j+1) = cosa;
    }
    // std::cout << "matB =\n"<< matB << std::endl;
    // std::cout << "rotation =\n"<< rotation << std::endl;
    // std::cout << "K =\n"<< matB.transpose() * matB << std::endl;
    MatrixDouble Mass(2*nphi,2*nphi);
    Mass = rotation * matB.transpose() * matB * rotation.transpose() * WJ * fDensity * fArea;

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
    MatrixDouble rotation(fDimension*nphi,fDimension*nphi);
    MatrixDouble matB(fDimension,fDimension*nphi);
    rotation.setZero();
    matB.setZero();
    double cosa = data.fAxes0(0,0) / data.fJacA0;
    double sina = data.fAxes0(1,0) / data.fJacA0;
    double check = sina*sina+cosa*cosa;
    for (int j = 0; j < nphi; j++){
        matB(0,fDimension*j) = data.fPhi[j];
        rotation(2*j  ,2*j  ) = cosa;
        rotation(2*j+1,2*j  ) = sina;
        rotation(2*j  ,2*j+1) = -sina;
        rotation(2*j+1,2*j+1) = cosa;
    }

    auto vel = data.fDSolDt;
    auto acel = data.fDSolDDt;
    auto disp = data.fSol;
    auto dispPrev = data.fSolPrev;

    // std::cout << "rotation =\n"<< rotation << std::endl;
    auto aux = rotation * matB.transpose() * WJ * fDensity * fArea;
    // std::cout << "Aux  = " << aux << std::endl;
    for (size_t i = 0; i < nphi; i++){
        auto aux = (disp/(fBeta * fTimeStep * fTimeStep) + 
                    vel/(fBeta*fTimeStep) + 
                    (1./(2.*fBeta)-1.) * acel) +
                    (disp*fGamma/(fBeta * fTimeStep) +
                    vel * (fGamma/fBeta - 1.) +
                    acel * fTimeStep * (fGamma/(2.*fBeta)-1.)) * fDamping;
        Rhs += rotation * matB.transpose() * aux * WJ * fDensity * fArea;
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
                auto dispPrev = cmesh->NodeVec()[inode]->PrevSolution();
                auto disp = cmesh->NodeVec()[inode]->Solution();
                VecDouble acelUpdated(2), velUpdated(2);
                acelUpdated = (disp-dispPrev)/(fBeta*fTimeStep*fTimeStep) -
                                velPrev/(fBeta*fTimeStep) -
                                (1./(2.*fBeta) - 1.) * acelPrev; 

                //Update Velocity
                velUpdated = velPrev + (1.-fGamma)*fTimeStep*acelPrev + fGamma*fTimeStep*acelUpdated;

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