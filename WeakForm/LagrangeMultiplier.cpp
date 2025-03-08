#include "LagrangeMultiplier.h"

LagrangeMultiplier::LagrangeMultiplier(int matid, int nstate, double scale) : WeakForm() {
    this->fMatId = matid;
    fNState = nstate;
    fScale = scale;
};

void LagrangeMultiplier::ComputeStiffness(int &index, IntPointData &leftdata, IntPointData &rightdata, MatrixDouble &Stiffness){
    if (!leftdata.fNeedsSol){
        leftdata.fNeedsSol = true;
        leftdata.fSol.resize(fNState);
    }
    if (!rightdata.fNeedsSol){
        rightdata.fNeedsSol = true;
        rightdata.fSol.resize(fNState);
    }
    double WJ = leftdata.fWeight * leftdata.fJacA0 * fScale;
    double nphi = leftdata.fPhi.size();
    
    Stiffness += WJ * leftdata.fPhi * rightdata.fPhi.transpose();
   
}

void LagrangeMultiplier::ComputeResidual(int &index, IntPointData &leftdata, IntPointData &rightdata, VecDouble &Rhs){

    double WJ = leftdata.fWeight * leftdata.fJacA0 * leftdata.fWeightFunction[index];
    double nphiL = leftdata.fPhi.size();
    double nphiR = rightdata.fPhi.size();
    VecDouble forcingF(1);
    VecDouble x_ = leftdata.fX;
    if (fForceFunction) fForceFunction(x_,forcingF);

    for (int i = 0; i < nphiL; i++){
        Rhs[i] -= WJ * leftdata.fPhi[i] * rightdata.fSol[0];
    }
    for (int i = 0; i < nphiR; i++){
        Rhs[nphiL+i] += (forcingF[0] * rightdata.fPhi[i] - rightdata.fPhi[i] * leftdata.fSol[0]) * WJ;
    }
};



int LagrangeMultiplier::VariableIndex(const std::string &name) const{
    
    if(!strcmp("LagrangeMultiplier",name.c_str()))    return 1;
    
    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int LagrangeMultiplier::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
        return fNState;

    default:
        PanicButton();
        return -1;
    }
};

void LagrangeMultiplier::Solution(IntPointData &data, int var, VecDouble &Sol) {

    //Solution
    if (var == 1){
        for (int i = 0; i < fNState; i++){
            Sol[i] = data.fSol[i];
        }
        return;
    };

}; 