#include "PhaseField.h"

//Class constructor
PhaseField::PhaseField(int matid, int dim, double kappa, double eta, double dTime) : WeakForm()
{
    fKappa = kappa;
    fEta = eta;
    fTimeStep = dTime;
    fIntegScheme = IEuler;
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
    double WJ = data.fWeight * data.fJacA0;
    int nphi = data.fPhi.size();

    // TODO: Use fWeightFunction instead of data.fJ.
    const auto sol = data.fSol(0, 0);

    const auto y = 1e-1*-fM * (data.fInterpWeightFunction + fGamma*(4*sol*sol - 6*sol + 2));//(1.0 - sol) * (sol - 0.5 - 30.0 * fEta * data.fJ * (1-sol)*sol);
    for (int i = nphi; i-- ; ){
        for (int j = nphi; j-- ; ){
            for (int k = fDimension; k--;  ){
                Stiffness(i,j) += fM*fGamma*fKsi*data.fDPhiX0(k,i) * data.fDPhiX0(k,j) * WJ;
            }
            Stiffness(i,j) += data.fPhi[i] * data.fPhi[j] * y * WJ;
        };
    };
}

void PhaseField::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){

    ComputeStiffnessStatic(index, data, Stiffness);

    const int nphi = data.fPhi.size();
    const double WJ = data.fWeight * data.fJacA0;
    MatrixDouble Mass(nphi,nphi);
    Mass.setZero();
    for (size_t i = 0; i < nphi; i++){
        for (size_t j = 0; j < nphi; j++){
            Mass(i,j) += data.fPhi[i] * data.fPhi[j] *  WJ;
        }
    }

    switch (fIntegScheme)
    {
    case IEuler:
        Stiffness += (1. / (fTimeStep)) * Mass;
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
    double WJ = data.fWeight * data.fJacA0;

    VecDouble forcingF(1);
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    const auto sol = data.fSol(0, 0);

    // TODO: Use fWeightFunction instead of data.fJ.
    const auto y = 1e-1*-fM * (data.fInterpWeightFunction + fGamma*(4*sol*sol - 6*sol + 2));//(1.0 - sol) * (sol - 0.5 - 30.0*fEta * data.fJ * (1-sol)*sol);

    for (int i = nphi; i--; ){
        double shapeFi = data.fPhi[i];

        //Matrix residual
        double K = 0.;
        for (int l=fDimension; l--; ) K += fM*fGamma*fKsi*data.fDPhiX0(l,i) * data.fDSolDx(0,l);

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
    const double WJ = data.fWeight * data.fJacA0;

    // int (R) -> int (w[dphi/dt + ...])
    // R_i+1 = R_i + grad(R) * delta_phi
    // dphi/dt -> int(w * dphi/dt)

    switch (fIntegScheme)
    {
    case IEuler:
        for (size_t i = 0; i < nphi; i++)
        {
            Rhs[i] -= ((data.fSol[0] - data.fSolPrev[0]) / fTimeStep) * data.fPhi[i] * WJ;
        }
        break;

    default:
        PanicButton();
        break;
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
    if(!strcmp("TimeDerivative",name.c_str()))  return 6;
    if(!strcmp("WeightFunction",name.c_str()))  return 7;
    
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
    case 6:
    case 7:
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

    //Velocity
    if (var == 6){
        Sol[0] = data.fDSolDt[0];
        return;
    };

    //WeightFunction
    if (var == 7){
        Sol[0] = data.fInterpWeightFunction;
        return;
    };
    
    
};

void PhaseField::UpdateTimeDerivatives(CompMesh *cmesh)
{
    switch (fIntegScheme)
    {
    case IEuler:
        for (int64_t inode = 0; inode < cmesh->NNodes(); inode++)
        {
            VecDouble dispPrev = cmesh->NodeVec()[inode]->PrevSolution();
            VecDouble disp = cmesh->NodeVec()[inode]->Solution();
            VecDouble velUpdated = (disp - dispPrev)/fTimeStep;

            cmesh->NodeVec()[inode]->SetDSolutionDTime(0, velUpdated[0]);
            cmesh->NodeVec()[inode]->SetPreviousSolution(0, disp[0]);
        }
        break;
    default:
        PanicButton();
        break;
    }
}
