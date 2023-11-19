
#include "PlasticityModel.h"


void PlasticityModel::BuildStressTensor(VecDouble &Sol, MatrixDouble &Stress){

    if (Sol.size() == 1){
        Stress(0,0) = Sol[0];
    } else if (Sol.size() == 3){
        Stress(0,0) = Sol[0];
        Stress(1,1) = Sol[1];
        Stress(0,1) = Stress(1,0) = Sol[2];
    } else if (Sol.size() == 6){
        Stress(0,0) = Sol[0];
        Stress(1,1) = Sol[1];
        Stress(2,2) = Sol[2];
        Stress(0,1) = Stress(1,0) = Sol[3];
        Stress(0,2) = Stress(2,0) = Sol[4];
        Stress(1,2) = Stress(2,1) = Sol[5];
    } else {
        PanicButton();
    }


}

void PlasticityModel::ComputePrincipalStress(MatrixDouble &Stress, VecDouble &PrincipalS){
    if(Stress.rows() == 1){
        PrincipalS[0] = Stress(0,0);
    } else {
        EigenSolver<MatrixDouble> solver(Stress,EigenvaluesOnly);
        PrincipalS = solver.eigenvalues().real();
    }
}