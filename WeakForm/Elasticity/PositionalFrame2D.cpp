#include "PositionalFrame2D.h"


PositionalFrame2D::PositionalFrame2D(int matid, double young, double inertia, double area, double height) : WeakForm() {
    this->fMatId = matid;
    fDimension = 1;
    fNState = 3;
    fYoungModulus = young;
    fInertia = inertia;
    fArea = area;
    fHeight = height;
};


void PositionalFrame2D::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    if (!data.fNeedsDSol || !data.fNeedsSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
        data.fNeedsSol = true;
        data.fSol.resize(fNState);
    }

    //Tangent vector
    VecDouble tangentVersor(2);
    double normTangent = sqrt(data.fAxes0(0,0)*data.fAxes0(0,0) + data.fAxes0(1,0)*data.fAxes0(1,0));
    tangentVersor[0] = -data.fAxes0(1,0)/normTangent;
    tangentVersor[1] =  data.fAxes0(0,0)/normTangent;
    double theta0 = atan2(tangentVersor[1],tangentVersor[0]);


    // Variable initialization
    double Length = 2.*data.fJacA0;
    double cosa = data.fAxes0(0,0) / data.fJacA0;
    double sina = data.fAxes0(1,0) / data.fJacA0;
    MatrixDouble rotation(6,6);
    rotation.setZero();

    //Local Stiffness
    Stiffness(0,0) = (fYoungModulus * fArea) / (Length);
    Stiffness(0,1) = Stiffness(1,0) = Stiffness(0,2) = Stiffness(2,0) = 0.;
    Stiffness(0,3) = Stiffness(3,0) = -(fYoungModulus * fArea) / (Length);
    Stiffness(0,4) = Stiffness(4,0) = Stiffness(0,5) = Stiffness(5,0) = 0.;

    Stiffness(1,1) = (12. * (fYoungModulus * fInertia) / (Length * Length * Length));
    Stiffness(1,2) = Stiffness(2,1) = (6. * (fYoungModulus * fInertia) / (Length * Length));
    Stiffness(1,3) = Stiffness(3,1) = 0.;
    Stiffness(1,4) = Stiffness(4,1) = -(12. * (fYoungModulus * fInertia) / (Length * Length * Length));
    Stiffness(1,5) = Stiffness(5,1) = (6. * (fYoungModulus * fInertia) / (Length * Length));

    Stiffness(2,2) = (4. * (fYoungModulus * fInertia) / Length);
    Stiffness(2,3) = Stiffness(3,2) = 0.;
    Stiffness(2,4) = Stiffness(4,2) = -(6. * (fYoungModulus * fInertia) / (Length * Length));
    Stiffness(2,5) = Stiffness(5,2) = (2. * (fYoungModulus * fInertia) / Length);

    Stiffness(3,3) = (fYoungModulus * fArea) / (Length);
    Stiffness(3,4) = Stiffness(4,3) = Stiffness(3,5) = Stiffness(5,3) = 0.;

    Stiffness(4,4) = (12. * (fYoungModulus * fInertia) / (Length * Length * Length));
    Stiffness(4,5) = Stiffness(5,4) = -(6. * (fYoungModulus * fInertia) / (Length * Length));

    Stiffness(5,5) = (4. * (fYoungModulus * fInertia) / Length);

    //Rotation matrix
    for (int j = 0; j < 2; j++){
        rotation(3*j  ,3*j  ) = cosa;
        rotation(3*j+1,3*j  ) = -sina;
        rotation(3*j  ,3*j+1) = sina;
        rotation(3*j+1,3*j+1) = cosa;
        rotation(3*j+2,3*j+2) = 1.;
    }

    Stiffness = rotation.transpose() * Stiffness * rotation;

}

void PositionalFrame2D::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    int nphi = 4;
    double WJ = data.fWeight * data.fJacA0;
    auto force = fForceFunction;
    VecDouble forcingF(3);
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    VecDouble phi,dphi,d2phi,d3phi;
    VecDouble Rhsaux(6);
    Rhsaux.setZero();

    HermiteFunction(data.fAdimCoord[0],phi,dphi,d2phi,d3phi);
    Rhsaux[0] = forcingF[0] * data.fPhi[0] * WJ;
    Rhsaux[3] = forcingF[0] * data.fPhi[1] * WJ;

    Rhsaux[1] = forcingF[1] * phi[0] * WJ;
    Rhsaux[4] = forcingF[1] * phi[2] * WJ;

    Rhsaux[2] = forcingF[1] * phi[1] * WJ * 0.5;
    Rhsaux[5] = forcingF[1] * phi[3] * WJ * 0.5;

    double cosa = data.fAxes0(0,0) / data.fJacA0;
    double sina = data.fAxes0(1,0) / data.fJacA0;
    MatrixDouble rotation(6,6);
    rotation.setZero();
    //Rotation matrix
    for (int j = 0; j < 2; j++){
        rotation(3*j  ,3*j  ) = cosa;
        rotation(3*j+1,3*j  ) = -sina;
        rotation(3*j  ,3*j+1) = sina;
        rotation(3*j+1,3*j+1) = cosa;
        rotation(3*j+2,3*j+2) = 1.;
    }

    Rhs += rotation.transpose() * Rhsaux;
    
};

void PositionalFrame2D::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();
}

int PositionalFrame2D::VariableIndex(const std::string &name) const{
    
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

int PositionalFrame2D::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
        return 3;
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

void PositionalFrame2D::Solution(IntPointData &data, int var, VecDouble &Sol) {

    auto force = fForceFunction;
    VecDouble forcingF(3);
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    VecDouble phi,dphi,d2phi,d3phi;

    HermiteFunction(data.fAdimCoord[0],phi,dphi,d2phi,d3phi);
    
    //Displacement
    if (var == 1){
        Sol[0] = data.fSol[0];
        Sol[1] = data.fSol[1];
        return;
    };

    //Rotation
    if (var == 2){
        Sol[0] = data.fSol[2];
        return;
    };

}; 

void PositionalFrame2D::HermiteFunction(double &ksi, VecDouble &phi, VecDouble &dphi, VecDouble &d2phi, VecDouble &d3phi){
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