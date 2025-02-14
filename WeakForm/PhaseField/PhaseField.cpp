#include "PhaseField.h"

//Class constructor
PhaseField::PhaseField(int matid, int dim, double dTime) : WeakForm()
{
    fTimeStep = dTime;
    fIntegScheme = EEuler;
    fMatId = matid;
    fDimension = dim;
    fNState = 1;
};

void PhaseField::ComputeStiffnessStatic(int &index, IntPointData &data, MatrixDouble &Stiffness)
{
    if (!data.fNeedsDSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(1,fDimension);
    }

    if (!data.fNeedsSol){
        data.fNeedsSol = true;
        data.fSol.resize(1);
    }

    // TODO: Remove fWeightFunction.
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    int nphi = data.fPhi.size();

    // TODO: Use fWeightFunction instead of data.fJ.
    const auto sol = data.fSol(0, 0);
    const auto y = (1.0 - sol) * (sol - 0.5 - 30.0 * fEta * data.fJ * (1-sol)*sol);
    for (int i = nphi; i-- ; ){
        for (int j = nphi; j-- ; ){
            for (int k = fDimension; k--;  ){
                Stiffness(i,j) += fKappa*data.fDPhiX0(k,i) * data.fDPhiX0(k,j) * WJ;
            }
            Stiffness(i,j) += data.fPhi[i] * data.fPhi[j] * y * WJ;
        };
    };
}

void PhaseField::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    ComputeStiffnessStatic(index, data, Stiffness);

    const int nphi = data.fPhi.size();
    const double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    MatrixDouble Mass(nphi,nphi);
    Mass.setZero();
    for (size_t i = 0; i < nphi; i++){
        for (size_t j = 0; j < nphi; j++){
            Mass(i,j) += data.fPhi[i] * data.fPhi[j] *  WJ;
        }
    }

    switch (fIntegScheme)
    {
    case EEuler:
        Stiffness += (1./(fTimeStep)) * Mass;    
        break;
    default:
        PanicButton();
        break;
    }
}

void PhaseField::ComputeResidualStatic(int &index, IntPointData &data, VecDouble &Rhs)
{
    auto force = fForceFunction;
    int nphi = data.fPhi.size();

    // TODO: Remove fWeightFunction.
    double WJ = data.fWeight * data.fJacA0  * data.fWeightFunction[index];

    VecDouble forcingF(1);
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    const auto sol = data.fSol(0, 0);

    // TODO: Use fWeightFunction instead of data.fJ.
    const auto y = (1.0 - sol) * (sol - 0.5 - 30.0*fEta * data.fJ * (1-sol)*sol);

    for (int i = nphi; i--; ){
        double shapeFi = data.fPhi[i];

        //Matrix residual
        double K = 0.;
        for (int l=fDimension; l--; ) K += fKappa*data.fDPhiX0(l,i) * data.fDSolDx(0,l);

        double phaseFieldRes = 0.0;

        phaseFieldRes = shapeFi * y * sol;

        //Source term
        double F = (forcingF[0]) * shapeFi;

        Rhs[i] += (-phaseFieldRes -K + F) * WJ;
    };
}

void PhaseField::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){
    ComputeResidualStatic(index, data, Rhs);
    
    const auto vel = data.fDSolDt;
    const int nphi = data.fPhi.size();
    const double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];

    // int (R) -> int (w[dphi/dt + ...])
    // R_i+1 = R_i + grad(R) * delta_phi
    // dphi/dt -> int(w * dphi/dt)
    for (size_t i = 0; i < nphi; i++){
        Rhs[i] += vel[0] * data.fPhi[i] * WJ;
    }
};

void PhaseField::ComputeError(IntPointData &data, VecDouble &errors){

    errors.resize(3);
    errors.setZero();

    VecDouble uExact(1);
    MatrixDouble DuExact(3,1);
    VecDouble x_ = data.fX;
    fExactSol(x_,uExact,DuExact);

    //Consider Arlequin weight function
    // uExact *= data.fWeightFunction[index];
    // DuExact *= data.fWeightFunction[index];
    // data.fSol *= data.fWeightFunction[index];
    // data.fDSolDx *= data.fWeightFunction[index];

    //L2 state variable
    errors[0] += (uExact[0]-data.fSol[0])*(uExact[0]-data.fSol[0]) * data.fWeight * data.fJacA0 ;
        
    //Semi H1 state variable
    for (int m = fDimension; m--; ){
        errors[1] += (DuExact(m,0)-data.fDSolDx(0,m))* (DuExact(m,0)-data.fDSolDx(0,m)) * data.fWeight * data.fJacA0;
    }

    //H1 state variable
    errors[2] = errors[0]+errors[1];
}


int PhaseField::VariableIndex(const std::string &name) const{
    
    if(!strcmp("Solution",name.c_str()))        return 1;
    if(!strcmp("Derivative",name.c_str()))      return 2;
    if(!strcmp("ExactSolution",name.c_str()))   return 3;
    if(!strcmp("ExactDerivative",name.c_str())) return 4;
    if(!strcmp("ForceFunction",name.c_str()))   return 5;
    
    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int PhaseField::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
    case 3:
    case 5:
        return 1;
    case 2:
    case 4:
        return 3;

    default:
        PanicButton();
        return -1;
    }
};

void PhaseField::Solution(IntPointData &data, int var, VecDouble &Sol) {

    //Solution
    if (var == 1){
        Sol[0] = data.fSol[0];
        return;
    };

    //Derivative
    if (var == 2){
        MatrixDouble aux = data.fAxes0 * data.fDSolDx.transpose();
        Sol[0] = aux(0,0);
        if (fDimension > 1) Sol[1] = aux(1,0);
        if (fDimension == 3) Sol[2] = aux(2,0);
        return;
    };

    VecDouble forcingF(1);
    VecDouble x_ = data.fX;
    if (fForceFunction) fForceFunction(x_,forcingF);

    VecDouble disp(3);
    MatrixDouble gradDisp(3,3);
    if (fExactSol) fExactSol(x_,disp,gradDisp);

    //Exact Solution
    if (var == 3){
        Sol[0] = disp[0];
        return;
    };
    
    //Exact Derivative
    if (var == 4){
        Sol[0] = gradDisp(0,0);
        Sol[1] = gradDisp(1,0);
        if (fDimension == 3) Sol[2] = gradDisp(2,0);
        return;
    };

    //ForceFunction
    if (var == 5){
        Sol[0] = forcingF[0];
        return;
    };
};

void PhaseField::UpdateTimeDerivatives(CompMesh *cmesh)
{
    switch (fIntegScheme)
    {
    case EEuler:
        for (int64_t inode = 0; inode < cmesh->NNodes(); inode++)
        {
            auto dispPrev = cmesh->NodeVec()[inode]->PrevSolution();
            auto disp = cmesh->NodeVec()[inode]->Solution();
            auto velUpdated = (disp - dispPrev)/fTimeStep;

            cmesh->NodeVec()[inode]->SetDSolutionDTime(0, velUpdated[0]);
        }
        break;
    default:
        PanicButton();
        break;
    }
}
