#include "MohrCoulomb.h"

void MohrCoulomb::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){

    Elasticity2D::ComputeStiffness(index,data,Stiffness);

    // double theta = ;


}

void MohrCoulomb::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    Elasticity2D::ComputeResidual(index,data,Rhs);

    
};

void MohrCoulomb::ComputeError(IntPointData &data, VecDouble &errors){
    
}

