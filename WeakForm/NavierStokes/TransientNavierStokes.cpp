#include "TransientNavierStokes.h"

TransientNavierStokes::TransientNavierStokes(int matid, int dim, double density, double viscosity) : NavierStokes(matid,dim,density,viscosity) {
    
};

void TransientNavierStokes::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    
    NavierStokes::ComputeStiffness(index,data,Stiffness);
    
}

void TransientNavierStokes::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    NavierStokes::ComputeResidual(index,data,Rhs);
    
};

void TransientNavierStokes::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet.\n";
    // PanicButton();
}

int TransientNavierStokes::VariableIndex(const std::string &name) const{
    PanicButton();
};

int TransientNavierStokes::NSolutionVariables(int var) const{
    PanicButton();
};

void TransientNavierStokes::Solution(IntPointData &data, int var, VecDouble &Sol){
    PanicButton();
};

void TransientNavierStokes::UpdateTimeDerivatives(CompMesh *cmesh){

};