#include "ElCoupling.h"


void ElCoupling::ComputeStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, MatrixDouble &Stiffness){

    int DIM = Mesh()->Dimension();
    //Velocity
    VecDouble u_(DIM);
    interpolateSolution(index, u_);

    //Velocity Derivatives
    MatrixDouble du_dx(DIM,DIM);
    interpolateSolDerivatives(dphi_dx, du_dx);

    //Lagrange Multiplier
    VecDouble lagM_(DIM);
    interpolateLagMultiplier(index, lagM_);

    //Lagrange Multiplier Derivatives
    MatrixDouble dL_dx(DIM,DIM);
    interpolateLagMultiplierDerivatives(dphi_dx, dL_dx);

    double k1 = fMeshVector[2]->getProblemParameters().getArlequinK1();
    double k2 = fMeshVector[2]->getProblemParameters().getArlequinK2();
        
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
            // // Lagrange multipliers residual
            // double L2 = lagM_[0] * Mesh()->getNumericalIntegration()-> phi_(i,index) * k1;

            // double H1 = 0.;
            // for (int l=DIM; l--; ) H1 += dphi_dx(i,l) * dL_dx(0,l) * k2;

            // lagrMultVector[i] += (L2 + H1) * WJ;

            // double L2u = una_[0] * Mesh()->getNumericalIntegration()-> phi_(i,index) * k1;

            // double H1u = 0.;
            // for (int l=DIM; l--; ) H1u += dphi_dx(i,l) * duna_dx(0,l) * k2;

            // rhsVector[i] += (L2u + H1u) * WJ;
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
                    };
                };
            };
            // Lagrange multipliers residual
            // for (int k = 0; k < DIM; k++){
            //     double L2 = lagM_[k] * fMesh->getNumericalIntegration()-> phi_(i,index) * k1;

            //     double H1 = 0.;
            //     for (int l=DIM; l--; ) H1 += dphi_dx(i,l) * dL_dx(k,l) * k2;
            //     for (int l=DIM; l--; ) H1 += dphi_dx(i,l) * dL_dx(l,k) * k2;

            //     lagrMultVector[DIM*i+k] += (L2 + H1) * WJ;

            //     double L2u = una_[k] * fMesh->getNumericalIntegration()-> phi_(i,index) * k1;

            //     double H1u = 0.;
            //     for (int l=DIM; l--; ) H1u += dphi_dx(i,l) * duna_dx(k,l) * k2;
            //     for (int l=DIM; l--; ) H1u += dphi_dx(i,l) * duna_dx(l,k) * k2;

            //     rhsVector[DIM*i+k] += (L2u + H1u) * WJ;
            // };
        };
    }






    // if (fMesh->ProbType() == ProblemType::EPoisson){
    //     for (int i = 0; i < fMesh->NElNodes(); i++){
    //         if (fMesh->NodeVec()[connect_[i]] -> getConstrains(0) == 1) {
    //             for (int j = 0; j < fMesh->NElNodes(); j++){
    //                 Stiffness(i,j) = 0.;
    //                 Stiffness(j,i) = 0.;
    //             };
    //             lagrMultVector[i] = 0.0;
    //         };
    //     };
    // } else {
    //     for (int i = 0; i < fMesh->NElNodes(); i++){
    //         for (int k = DIM; k--; ){
    //             if (fMesh->NodeVec()[connect_[i]] -> getConstrains(k) == 1) {
    //                 for (int j = 0; j < fMesh->NElNodes()*DIM; j++){
    //                     Stiffness(DIM*i+k,j) = 0.;
    //                     Stiffness(j,DIM*i+k) = 0.;
    //                 };
    //                 lagrMultVector[DIM*i+k] = 0.0;
    //             };
    //         };
    //     };
    // }
    



}


