#include "LinearBeam.h"


LinearBeam::LinearBeam(int matid, REAL young, REAL inertia) : WeakForm() {
    this->fMatId = matid;
    fDimension = 2;
    fNState = 2;
    fYoungModulus = young;
    fInertia = inertia;
    this->fType = WeakFormType::kLinearBeam;
};


void LinearBeam::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    if (!data.fNeedsDSol || !data.fNeedsSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
        data.fNeedsSol = true;
        data.fSol.resize(fNState);
    }

    REAL Length = 2.*data.fJacA0;

    //Local Stiffness
    Stiffness(0,0) = (12 * (fYoungModulus * fInertia) / (Length * Length * Length));
    Stiffness(0,1) = Stiffness(1,0) = (6 * (fYoungModulus * fInertia) / (Length * Length));
    Stiffness(0,2) = Stiffness(2,0) = -(12 * (fYoungModulus * fInertia) / (Length * Length * Length));
    Stiffness(0,3) = Stiffness(3,0) = (6 * (fYoungModulus * fInertia) / (Length * Length));

    Stiffness(1,1) = (4 * (fYoungModulus * fInertia) / Length);
    Stiffness(1,2) = Stiffness(2,1) = -(6 * (fYoungModulus * fInertia) / (Length * Length));
    Stiffness(1,3) = Stiffness(3,1) = (2 * (fYoungModulus * fInertia) / Length);    

    Stiffness(2,2) = (12 * (fYoungModulus * fInertia) / (Length * Length * Length));
    Stiffness(2,3) = Stiffness(3,2) = -(6 * (fYoungModulus * fInertia) / (Length * Length));

    Stiffness(3,3) = (4 * (fYoungModulus * fInertia) / Length);

}

void LinearBeam::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    int nphi = 4;
    REAL WJ = data.fWeight * data.fJacA0;
    auto force = fForceFunction;
    VecDouble forcingF(fDimension);
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    VecDouble phi,dphi,d2phi,d3phi;

    HermiteFunction(data.fAdimCoord[0],phi,dphi,d2phi,d3phi);
    
    for (int i = 4; i--; ){
        //External force
        REAL corr = i%2 == 0 ? 1. : data.fJacA0; // in this case the jacobian needs to be squared
        REAL Fx = forcingF[0] * phi[i] * corr;
        Rhs[i] += Fx * WJ;
    };
    
};

void LinearBeam::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();
}

int LinearBeam::VariableIndex(const std::string &name) const{
    
    if(!strcmp("Displacement",name.c_str()))           return 1;
    if(!strcmp("Rotation",name.c_str()))               return 2;
    if(!strcmp("ExactDisplacement",name.c_str()))      return 3;
    if(!strcmp("ExactRotation",name.c_str()))          return 4;
    if(!strcmp("ShearForce",name.c_str()))             return 5;
    if(!strcmp("BendingMoment",name.c_str()))          return 6;
    if(!strcmp("Load",name.c_str()))                   return 7;

    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int LinearBeam::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
    case 3:
    case 5:
    case 2:
    case 4:
    case 6:
    case 7:
        return 1;

    default:
        PanicButton();
        return -1;
    }
};

void LinearBeam::Solution(IntPointData &data, int var, VecDouble &Sol) {

    auto force = fForceFunction;
    VecDouble forcingF(fDimension);
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    VecDouble phi,dphi,d2phi,d3phi;

    HermiteFunction(data.fAdimCoord[0],phi,dphi,d2phi,d3phi);
    
    //Displacement
    if (var == 1){
        Sol[0] = data.fSol[0];
        return;
    };

    //Rotation
    if (var == 2){
        Sol[0] = data.fSol[1];
        return;
    };

    //ShearForce
    if (var == 5){
        Sol[0] = data.fSol[1];
        return;
    };

    //BendingMoment
    if (var == 6){
        if(data.fAdimCoord[0] == -1.){
            Sol[0] = data.fSol[1]*d2phi[1];
        } else{
            Sol[0] = data.fSol[1]*d2phi[3];
        }
        return;
    };


}; 


void LinearBeam::HermiteFunction(REAL &ksi, VecDouble &phi, VecDouble &dphi, VecDouble &d2phi, VecDouble &d3phi){
    phi.resize(4);
    phi[0] = 0.5 - 0.75*ksi + 0.25*ksi*ksi*ksi;
    phi[1] = 0.25 - 0.25*ksi - 0.25*ksi*ksi + 0.25*ksi*ksi*ksi;
    phi[2] = 0.5 + 0.75*ksi - 0.25*ksi*ksi*ksi;
    phi[3] = -0.25 - 0.25*ksi + 0.25*ksi*ksi + 0.25*ksi*ksi*ksi;
    dphi.resize(4);
    dphi[0] = -0.75 + 0.75*ksi*ksi;
    dphi[1] = -0.25 - 0.5*ksi + 0.75*ksi*ksi;
    dphi[2] = 0.75 - 0.75*ksi*ksi;
    dphi[3] = -0.25 + 0.5*ksi + 0.75*ksi*ksi;
    d2phi.resize(4);
    d2phi[0] = 1.5*ksi;
    d2phi[1] = -0.5 + 1.5*ksi;
    d2phi[2] = -1.5*ksi;
    d2phi[3] = 0.5 + 1.5*ksi;
    d3phi.resize(4);
    d3phi[0] = 1.5;
    d3phi[1] = 1.5;
    d3phi[2] = -1.5;
    d3phi[3] = 1.5;

}