
#include "PlasticityModel.h"

PlasticityModel::PlasticityModel(WeakForm *elast){
    fElasticModel = elast;
    fNState = elast->NState();
    fMatId = elast->Id();
    int var = fElasticModel->VariableIndex("Stress");
    int nsol = fElasticModel->NSolutionVariables(var);
    fRealDimension = 0;
    if (nsol == 1) fRealDimension = 1;
    if (nsol == 3) fRealDimension = 2;
    if (nsol == 6) fRealDimension = 3;

    fDimension = elast->Dimension();
};

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

void PlasticityModel::TensorToVoigt(MatrixDouble &tensor, VecDouble &voigt){

    switch (fRealDimension)
    {
    case 1:
        if (voigt.size() != 1) voigt.resize(1);
        voigt[0] = tensor(0,0);
        break;
    case 2:
        if (voigt.size() != 3) voigt.resize(3);
        voigt[0] = tensor(0,0);
        voigt[1] = tensor(1,1);
        voigt[2] = tensor(0,1);
        break;
    case 3:
        if (voigt.size() != 6) voigt.resize(6);
        voigt[0] = tensor(0,0);
        voigt[1] = tensor(1,1);
        voigt[2] = tensor(2,2);
        voigt[3] = tensor(0,1);
        voigt[4] = tensor(0,2);
        voigt[5] = tensor(1,2);
        break;
    
    default:
        PanicButton();
        break;
    }
};

void PlasticityModel::VoigtToTensor(MatrixDouble &tensor, VecDouble &voigt){
    switch (fRealDimension)
    {
    case 1:
        if (tensor.rows() != 1) tensor.resize(1,1);
        tensor(0,0) = voigt[0];
        break;
    case 2:
        if (tensor.rows() != 2) tensor.resize(2,2);
        tensor(0,0) = voigt[0];
        tensor(1,1) = voigt[1];
        tensor(0,1) = tensor(1,0) = voigt[2];
        break;
    case 3:
        if (tensor.rows() != 3) tensor.resize(3,3);
        tensor(0,0) = voigt[0];
        tensor(1,2) = voigt[1];
        tensor(2,2) = voigt[2];
        tensor(0,1) = tensor(1,0) = voigt[3];
        tensor(0,2) = tensor(2,0) = voigt[4];
        tensor(1,2) = tensor(2,1) = voigt[5];
        break;

    default:
        PanicButton();
        break;
    }

};