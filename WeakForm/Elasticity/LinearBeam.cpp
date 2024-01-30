#include "LinearBeam.h"


LinearBeam::LinearBeam(int matid, int dim, double young, double inertia) : WeakForm() {
    this->fMatId = matid;
    fDimension = dim;
    fNState = dim;
    fYoungModulus = young;
    fInertia = inertia;
};


void LinearBeam::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    // if (!data.fNeedsDSol || !data.fNeedsSol){
    //     data.fNeedsDSol = true;
    //     data.fDSolDx.resize(fNState,fDimension);
    //     data.fNeedsSol = true;
    //     data.fSol.resize(fNState);
    // }
    
    // double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    // int nphi = data.fPhi.size();
    // double elementLenght = 2.*data.fJacA0;
    // double K = fYoungModulus * fArea / elementLenght;


    // MatrixDouble rotation(fDimension*nphi,fDimension*nphi);
    // MatrixDouble matB(fDimension,fDimension*nphi);
    // rotation.setZero();
    // matB.setZero();
    // double cosa = data.fAxes0(0,0) / data.fJacA0;
    // double sina = data.fAxes0(1,0) / data.fJacA0;
    // double check = sina*sina+cosa*cosa;
    // for (int j = 0; j < nphi; j++){
    //     // for (int i = 0; i < fDimension; i++){
    //         matB(0,fDimension*j) = data.fDPhiX0(j,0);
    //     // }
    //     rotation(2*j  ,2*j  ) = cosa;
    //     rotation(2*j+1,2*j  ) = sina;
    //     rotation(2*j  ,2*j+1) = -sina;
    //     rotation(2*j+1,2*j+1) = cosa;
    // }
    // // std::cout << "matB =\n"<< matB << std::endl;
    // // std::cout << "rotation =\n"<< rotation << std::endl;
    // // std::cout << "K =\n"<< matB.transpose() * matB << std::endl;
    // Stiffness += rotation * matB.transpose() * matB * rotation.transpose() * WJ * elementLenght * K;

}

void LinearBeam::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    // int nphi = data.fPhi.size();
    // double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    // auto force = fForceFunction;
    // VecDouble forcingF(fDimension);
    // VecDouble x_ = data.fX;
    // if (force) force(x_,forcingF);
    
    // double elementLenght = 2.*data.fJacA0;
    // double K = fYoungModulus * fArea / elementLenght;

    // MatrixDouble rotation(fDimension*nphi,fDimension*nphi);
    // MatrixDouble matB(fDimension,fDimension*nphi);
    // rotation.setZero();
    // matB.setZero();
    // double cosa = data.fAxes0(0,0) / data.fJacA0;
    // double sina = data.fAxes0(1,0) / data.fJacA0;
    // double check = sina*sina+cosa*cosa;
    // for (int j = 0; j < nphi; j++){
    //     // for (int i = 0; i < fDimension; i++){
    //         matB(0,fDimension*j) = data.fDPhiX0(j,0);
    //     // }
    //     rotation(2*j  ,2*j  ) = cosa;
    //     rotation(2*j+1,2*j  ) = sina;
    //     rotation(2*j  ,2*j+1) = -sina;
    //     rotation(2*j+1,2*j+1) = cosa;
    // }

    // VecDouble sol(2);// = data.fSol;
    // sol[0] = -data.fDSolDx(1,0)*sina;
    // sol[1] = +data.fDSolDx(0,0)*cosa;

    // // std::cout << "rotation =\n"<< rotation << std::endl;
    // Rhs += rotation * matB.transpose() * sol * WJ * elementLenght * K;

    // for (int i = nphi; i--; ){
    //     double shapeFi = data.fPhi[i];
    //     //External force
    //     double Fx = forcingF[0] * shapeFi;
    //     double Fy = forcingF[1] * shapeFi;
    //     Rhs[2*i  ] += Fx * WJ;
    //     Rhs[2*i+1] += Fy * WJ;
    // };
    // // std::cout << "RHS = " << Rhs << std::endl;
    
};

void LinearBeam::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();
}

int LinearBeam::VariableIndex(const std::string &name) const{
    
    if(!strcmp("Displacement",name.c_str()))           return 1;
    if(!strcmp("Rotation",name.c_str()))               return 2;
    if(!strcmp("ExactDisplacement",name.c_str()))      return 3;
    if(!strcmp("Rotation",name.c_str()))               return 4;
    if(!strcmp("ShearForce",name.c_str()))             return 5;
    if(!strcmp("BendingMoment",name.c_str()))          return 6;

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
        return 3;
    case 2:
    case 4:
        return 1;

    default:
        PanicButton();
        return -1;
    }
};

void LinearBeam::Solution(IntPointData &data, int var, VecDouble &Sol) {

    // //Displacement
    // if (var == 1){
    //     Sol[0] = data.fSol[0];
    //     Sol[1] = data.fSol[1];
    //     Sol[2] = 0.;
    //     return;
    // };

    // //NormalStress
    // if (var == 2){
    //     double cosa = data.fAxes0(0,0) / data.fJacA0;
    //     double sina = data.fAxes0(1,0) / data.fJacA0;
    //     Sol[0] = fYoungModulus * (data.fDSolDx(1,0)*sina - data.fDSolDx(0,0)*cosa) ;
    //     return;
    // };

    // VecDouble forcingF(fDimension);
    // VecDouble x_ = data.fX;
    // if (fForceFunction) fForceFunction(x_,forcingF);

    // VecDouble disp(fDimension);
    // MatrixDouble gradDisp(fDimension,3);
    // if (fExactSol) fExactSol(x_,disp,gradDisp);

    // //Exact Displacement
    // if (var == 3){
    //     Sol[0] = disp[0];
    //     Sol[1] = disp[1];
    //     Sol[2] = 0.;
    //     return;
    // };
    
    // //Exact Sigma X
    // if (var == 4){
    //     double epsilon = gradDisp.norm();
    //     Sol[0] = fYoungModulus * epsilon;
    //     return;
    // };

    // //Exact Force
    // if (var == 5){
    //     Sol[0] = forcingF[0];
    //     Sol[1] = forcingF[1];
    //     Sol[2] = 0.;
    //     return;
    // };

}; 