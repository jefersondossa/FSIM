#include "PositionalTruss.h"

PositionalTruss::PositionalTruss(int matid, int dim, REAL young, REAL area) : WeakForm() {
    this->fMatId = matid;
    fDimension = dim;
    fNState = dim;
    fYoungModulus = young;
    fArea = area;
    this->fType = WeakFormType::kPositionalTruss;
};

void PositionalTruss::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    if (!data.fNeedsDSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
    }

    int nphi = data.fPhi.size();
    REAL WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    
    //Green-Lagrange strain tensor
    REAL initLenght = 2.*data.fJacA0;
    REAL currLenght = 2.*data.fJacA1;
    REAL E = 0.5 * (currLenght*currLenght/(initLenght*initLenght) - 1.);
    
    //Piola-Kirchhoff stress tensor
    REAL S = fYoungModulus * E;

    VecDouble DeltaY = (data.fAxes1 * 2.) / initLenght;

    // stiffness matrix
    for (int a = 0; a < nphi; a++){
        for (int k = 0; k < fDimension; k++){
            //element tangent matrix
            for (int b = 0; b < nphi; b++){
                for (int l = 0; l < fDimension; l++){
                    REAL dkronecker = 0.;
                    if (k==l) dkronecker = 1.;
                    REAL aux = pow(-1.,a+1) * pow(-1.,b+1) * (fYoungModulus*DeltaY[k]*DeltaY[l] + S*dkronecker) * (fArea / initLenght);

                    Stiffness(fDimension*a+k,fDimension*b+l) += aux * data.fWeight / 2.;
                }
            }
        }
    }


}

void PositionalTruss::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    int nphi = data.fPhi.size();
    
    //Green-Lagrange strain tensor
    REAL initLenght = 2.*data.fJacA0;
    REAL currLenght = 2.*data.fJacA1;
    REAL E = 0.5 * (currLenght*currLenght/(initLenght*initLenght) - 1.);
    
    //Piola-Kirchhoff stress tensor
    REAL S = fYoungModulus * E;

    VecDouble DeltaY = (data.fAxes1 * 2.);

    // Residual vector
    for (int b = 0; b < nphi; b++){
        for (int k = 0; k < fDimension; k++){
            //Internal force
            REAL aux = fArea * S * pow(-1.,b+1) * DeltaY[k] / currLenght;
            
            Rhs[fDimension*b + k] -= aux * data.fWeight / 2.;
            
        }
    }
};

void PositionalTruss::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();

}


int PositionalTruss::VariableIndex(const std::string &name) const{
    
    if(!strcmp("Displacement",name.c_str()))           return 1;
    if(!strcmp("Stress",name.c_str()))           return 2;
    if(!strcmp("ExactDisplacement",name.c_str()))      return 3;
    if(!strcmp("ExactStress",name.c_str()))      return 4;
    if(!strcmp("ExactForce",name.c_str()))             return 5;

    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int PositionalTruss::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
    case 3:
    case 5:
        return 3;
    case 2:
    case 4:
        return 1;

    default:
        PanicButton();
        return -1;
    }
};

void PositionalTruss::Solution(IntPointData &data, int var, VecDouble &Sol) {

    //Displacement
    if (var == 1){
        Sol[0] = data.fSol[0];
        Sol[1] = data.fSol[1];
        Sol[2] = 0.;
        return;
    };

    //NormalStress
    if (var == 2){
        REAL cosa = data.fAxes0(0,0) / data.fJacA0;
        REAL sina = data.fAxes0(1,0) / data.fJacA0;
        Sol[0] = fYoungModulus * (data.fDSolDx(1,0)*sina - data.fDSolDx(0,0)*cosa) ;
        return;
    };

    VecDouble forcingF(fDimension);
    VecDouble x_ = data.fX;
    if (fForceFunction) fForceFunction(x_,forcingF);

    VecDouble disp(fDimension);
    MatrixDouble gradDisp(fDimension,3);
    if (fExactSol) fExactSol(x_,disp,gradDisp);

    //Exact Displacement
    if (var == 3){
        Sol[0] = disp[0];
        Sol[1] = disp[1];
        Sol[2] = 0.;
        return;
    };
    
    //Exact Sigma X
    if (var == 4){
        REAL epsilon = gradDisp.norm();
        Sol[0] = fYoungModulus * epsilon;
        return;
    };

    //Exact Force
    if (var == 5){
        Sol[0] = forcingF[0];
        Sol[1] = forcingF[1];
        Sol[2] = 0.;
        return;
    };

}; 

MatrixDouble &PositionalTruss::ConstitutiveMatrix(){
    MatrixDouble constitutive(1,1);
    constitutive(0,0) = fYoungModulus;
    return constitutive;
}