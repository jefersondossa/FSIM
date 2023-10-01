#include "ElCouplingLocal.h"


void ElCouplingLocal::ComputeStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, MatrixDouble &Stiffness){

    int DIM = Mesh()->Dimension();
    int DEG = Mesh()->GetDefaultOrder();

    double k1 = fMeshVector[1]->getProblemParameters().getArlequinK1();
    double k2 = fMeshVector[1]->getProblemParameters().getArlequinK2();
        
    double WJ = weight_ * djac_;

    if (Mesh()->ProbType() == ProblemType::EPoisson){
        for (int i = 0; i < Mesh()->NElNodes(); i++){
            for (int j = 0; j < Mesh()->NElNodes(); j++){
                double l2 = Mesh()->getNumericalIntegration()-> phi_(i,index) * Mesh()->getNumericalIntegration()-> phi_(j,index) * WJ * k1;
                // L2 COUPLING OPERATOR
                Stiffness(i,j) -= l2;
                for (int l = 0; l < DIM; l++){
                    //H1 COUPLING OPERATOR
                    double K = dphi_dx(i,l) * dphi_dx(j,l);
                    Stiffness(i,j) -= K * WJ * k2;
                };
            };
        };
    } else {
        for (int i = 0; i < Mesh()->NElNodes(); i++){
            for (int j = 0; j < Mesh()->NElNodes(); j++){
                double l2 = Mesh()->getNumericalIntegration()-> phi_(i,index) * Mesh()->getNumericalIntegration()-> phi_(j,index) * WJ * k1;
                for (int k = 0; k < DIM; k++){
                    // L2 COUPLING OPERATOR
                    Stiffness(DIM*i+k,DIM*j+k) -= l2;
                    for (int l = 0; l < DIM; l++){
                        //H1 COUPLING OPERATOR
                        double K = dphi_dx(i,l) * dphi_dx(j,k);
                        if (k==l) for (int m = DIM; m--; ) K += dphi_dx(i,m) * dphi_dx(j,m);

                        Stiffness(DIM*i+k,DIM*j+l) -= K * WJ * k2;
                    };//l
                };//k
            };//j
        };//i
    };//if Poisson
};


void ElCouplingLocal::ComputeResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, VecDouble &Rhs){

    int DIM = Mesh()->Dimension();
    int DEG = Mesh()->GetDefaultOrder();

    auto ellocal = fMeshVector[1]->ElementVec()[fLocalIndex];

    //Velocity
    VecDouble u_(DIM);
    ellocal->interpolateSolution(index, u_);

    //Velocity Derivatives
    MatrixDouble du_dx(DIM,DIM);
    ellocal->interpolateSolDerivatives(dphi_dx, du_dx);

    //The lagrange multiplier is the solution of meshvector[2]
    //Lagrange Multiplier
    VecDouble lagM_(DIM);
    interpolateSolution(index, lagM_);

    //Lagrange Multiplier Derivatives
    MatrixDouble dL_dx(DIM,DIM);
    interpolateSolDerivatives(dphi_dx, dL_dx);

    double k1 = fMeshVector[1]->getProblemParameters().getArlequinK1();
    double k2 = fMeshVector[1]->getProblemParameters().getArlequinK2();
        
    double WJ = weight_ * djac_;

    if (Mesh()->ProbType() == ProblemType::EPoisson){
        for (int i = 0; i < Mesh()->NElNodes(); i++){
            //Solution residual
            double L2u = u_[0] * Mesh()->getNumericalIntegration()-> phi_(i,index) * k1;
            double H1u = 0.;
            for (int l=DIM; l--; ) H1u += dphi_dx(i,l) * du_dx(0,l) * k2;
            Rhs[i] += (L2u + H1u) * WJ;
            
            // Lagrange multipliers residual
            double L2 = lagM_[0] * Mesh()->getNumericalIntegration()-> phi_(i,index) * k1;
            double H1 = 0.;
            for (int l=DIM; l--; ) H1 += dphi_dx(i,l) * dL_dx(0,l) * k2;
            Rhs[Mesh()->NLocDOF()+i] += (L2 + H1) * WJ;
        };
    } else {
        for (int i = 0; i < Mesh()->NElNodes(); i++){
            for (int k = 0; k < DIM; k++){
                //Solution residual
                double L2u = -u_[k] * Mesh()->getNumericalIntegration()-> phi_(i,index) * k1;
                double H1u = 0.;
                for (int l=DIM; l--; ) H1u += dphi_dx(i,l) * du_dx(k,l) * k2;
                for (int l=DIM; l--; ) H1u += dphi_dx(i,l) * du_dx(l,k) * k2;
                Rhs[DIM*i+k] += (L2u + H1u) * WJ;
                
                //Lagrange multipliers residual
                double L2 = lagM_[k] * Mesh()->getNumericalIntegration()-> phi_(i,index) * k1;
                double H1 = 0.;
                for (int l=DIM; l--; ) H1 += dphi_dx(i,l) * dL_dx(k,l) * k2;
                for (int l=DIM; l--; ) H1 += dphi_dx(i,l) * dL_dx(l,k) * k2;
                Rhs[Mesh()->NLocDOF()+DIM*i+k] += (L2 + H1) * WJ;
            };
        };
    }
}


void ElCouplingLocal::ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs){

    int DIM = Mesh()->Dimension();

    if (Mesh()->ProbType() == ProblemType::EPoisson){
        for (int i = 0; i < Mesh()->NElNodes(); i++){
            auto connectlocal = fMeshVector[1]->ElementVec()[fLocalIndex]->getConnectivity();
            if (fMeshVector[1]->NodeVec()[connectlocal[i]] -> getConstrains(0) == 1) {
                for (int j = 0; j < Mesh()->NElNodes(); j++){
                    Stiffness(i,j) = 0.;
                    // Stiffness(j,i) = 0.;
                };
                Rhs[i] = 0.0;
                // Rhs[Mesh()->NLocDOF() + i] = 0.0;
            };
        };
    } else {
        for (int i = 0; i < Mesh()->NElNodes(); i++){
            auto connectlocal = fMeshVector[1]->ElementVec()[fLocalIndex]->getConnectivity();
            for (int k = DIM; k--; ){
                if (fMeshVector[1]->NodeVec()[connectlocal[i]] -> getConstrains(k) == 1) {
                    for (int j = 0; j < Mesh()->NElNodes()*DIM; j++){
                        Stiffness(DIM*i+k,j) = 0.;
                        Stiffness(j,DIM*i+k) = 0.;
                    };
                    Rhs[DIM*i+k] = 0.0;
                    Rhs[Mesh()->NLocDOF() + DIM*i+k] = 0.0;
                };
            };
        };
    }




}

