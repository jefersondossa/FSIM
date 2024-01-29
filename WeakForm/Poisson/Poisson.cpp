#include "Poisson.h"

//Class constructor
Poisson::Poisson(int matid, int dim, int nState) : WeakForm() {
    fMatId = matid;
    fDimension = dim;
    fNState = nState;
};

void Poisson::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    if (!data.fNeedsDSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
    }

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    int nphi = data.fPhi.size();

    for (int i = nphi; i-- ; ){
        for (int j = nphi; j-- ; ){            
            for (int k = fDimension; k--;  ){
                Stiffness(i,j) += data.fDPhiX0(i,k) * data.fDPhiX0(j,k) * WJ;
            }
        };
    };
}

void Poisson::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    auto force = fForceFunction;
    int nphi = data.fPhi.size();

    double WJ = data.fWeight * data.fJacA0  * data.fWeightFunction[index];

    VecDouble forcingF(1);
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    for (int i = nphi; i--; ){
        double shapeFi = data.fPhi[i];

        //Matrix residual
        double K = 0.;
        for (int l=fDimension; l--; ) K += data.fDPhiX0(i,l) * data.fDSolDx(0,l);

        //Source term
        double F = (forcingF[0]) * shapeFi;

        Rhs[i] += (-K + F) * WJ;
    };
};

void Poisson::ComputeError(IntPointData &data, VecDouble &errors){

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


int Poisson::VariableIndex(const std::string &name) const{
    
    if(!strcmp("Solution",name.c_str()))        return 1;
    if(!strcmp("Derivative",name.c_str()))      return 2;
    if(!strcmp("ExactSolution",name.c_str()))   return 3;
    if(!strcmp("ExactDerivative",name.c_str())) return 4;
    if(!strcmp("ForceFunction",name.c_str()))   return 5;
    
    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int Poisson::NSolutionVariables(int var) const{
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

void Poisson::Solution(IntPointData &data, int var, VecDouble &Sol) {

    //Solution
    if (var == 1){
        Sol[0] = data.fSol[0];
        return;
    };

    //Derivative
    if (var == 2){
        Sol[0] = data.fDSolDx(0,0);
        if (fDimension > 1) Sol[1] = data.fDSolDx(0,1);
        if (fDimension == 3) Sol[2] = data.fDSolDx(0,2);
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