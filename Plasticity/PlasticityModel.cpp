
#include "PlasticityModel.h"
#include "ElasticTruss.h"
#include "Elasticity2D.h"
#include "PositionalTruss.h"
#include "ElasticityPositional2D.h"

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
    fNStressComponents = nsol;
    fDimension = elast->Dimension();

    ElasticTruss *truss = dynamic_cast<ElasticTruss* >(fElasticModel);
    Elasticity2D *mat2d = dynamic_cast<Elasticity2D* >(fElasticModel);
    PositionalTruss *postruss = dynamic_cast<PositionalTruss* >(fElasticModel);
    if (truss) {
        fYoungModulus = truss->YoungModulus();
        fPoissonRatio = 0.;
    }
    if (mat2d) {
        fYoungModulus = mat2d->YoungModulus();
        fPoissonRatio = mat2d->PoissonRatio();
        fPlaneStress = mat2d->PlaneState();
    }
    if (postruss) {
        fYoungModulus = postruss->YoungModulus();
        fPoissonRatio = 0.;
    }

    Tensor Id2;
    Id2.Identity();
    fId2xId2 = Id2.TensorProduct(Id2);
    
    std::vector<std::vector<std::vector<std::vector<double>>>> fIdSymmetric;
    fIdSymmetric.resize(3);
    for (int i = 0; i < 3; i++){
        fIdSymmetric[i].resize(3);
        for (int j = 0; j < 3; j++){
            fIdSymmetric[i][j].resize(3);
            for (int k = 0; k < 3; k++){
                fIdSymmetric[i][j][k].resize(3);
            }
        }
    }

    double dik,djl,dil,djk;

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                for (int l = 0; l < 3; l++){
                    dik = i == k ? 1. : 0.; 
                    djl = j == l ? 1. : 0.; 
                    dil = i == l ? 1. : 0.; 
                    djk = j == k ? 1. : 0.; 
                    fIdSymmetric[i][j][k][l] = 0.5 * (dik*djl + dil*djk);
                }
            }
        }
    }

    //4th order identity tensor
    MatrixDouble fIdentity4S(6,6);
    fIdentity4S.setZero();
    
    //Transform 4th order tensor into matrix: https://wiki.seg.org/wiki/Voigt_notation
    //Diagonal
    fIdentity4S(0,0) = fIdSymmetric[0][0][0][0];
    fIdentity4S(1,1) = fIdSymmetric[1][1][1][1];
    fIdentity4S(2,2) = fIdSymmetric[2][2][2][2];
    fIdentity4S(3,3) = fIdSymmetric[1][2][1][2];
    fIdentity4S(4,4) = fIdSymmetric[0][2][0][2];
    fIdentity4S(5,5) = fIdSymmetric[0][1][0][1];
    //1st line
    fIdentity4S(0,1) = fIdentity4S(1,0) = fIdSymmetric[0][0][1][1];
    fIdentity4S(0,2) = fIdentity4S(2,0) = fIdSymmetric[0][0][2][2];
    fIdentity4S(0,3) = fIdentity4S(3,0) = fIdSymmetric[0][0][1][2];
    fIdentity4S(0,4) = fIdentity4S(4,0) = fIdSymmetric[0][0][0][2];
    fIdentity4S(0,5) = fIdentity4S(5,0) = fIdSymmetric[0][0][0][1];
    //2nd line
    fIdentity4S(1,2) = fIdentity4S(2,1) = fIdSymmetric[1][1][2][2];
    fIdentity4S(1,3) = fIdentity4S(3,1) = fIdSymmetric[1][1][1][2];
    fIdentity4S(1,4) = fIdentity4S(4,1) = fIdSymmetric[1][1][0][2];
    fIdentity4S(1,5) = fIdentity4S(5,1) = fIdSymmetric[1][1][0][1];
    //3rd line
    fIdentity4S(2,3) = fIdentity4S(3,2) = fIdSymmetric[2][2][1][2];
    fIdentity4S(2,4) = fIdentity4S(4,2) = fIdSymmetric[2][2][0][2];
    fIdentity4S(2,5) = fIdentity4S(5,2) = fIdSymmetric[2][2][0][1];
    //4th line
    fIdentity4S(3,4) = fIdentity4S(4,3) = fIdSymmetric[1][2][0][2];
    fIdentity4S(3,5) = fIdentity4S(5,3) = fIdSymmetric[1][2][0][1];
    //5th line
    fIdentity4S(4,5) = fIdentity4S(5,4) = fIdSymmetric[0][2][0][1];

    fIdentity4Dev = fIdentity4S - fId2xId2/3.;

    if (fRealDimension == 2){
        if (fPlaneStress){
            fBulkModulus = fYoungModulus / (2. * (1.-fPoissonRatio));
        } else {
            fBulkModulus = fYoungModulus / (3. * (1.-2.*fPoissonRatio));
        }
        fShearModulus = fYoungModulus / (2. * (1.+fPoissonRatio));
    } else if (fRealDimension == 3) {
        fBulkModulus = fYoungModulus / (3. * (1.-2.*fPoissonRatio));
        fShearModulus = fYoungModulus / (2. * (1.+fPoissonRatio));
        
        std::cout << "please define the tensors for the desired dimension\n";
        PanicButton();
    }

    
};

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

void PlasticityModel::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs, Tensor &Stress){

    // fElasticModel->ComputeResidual(index,data,Rhs);
    // return;
    int nphi = data.fPhi.size();

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    MatrixDouble matB(3,2*nphi);
    matB.setZero();
    for (int j = 0; j < nphi; j++){
        matB(0,2*j  ) = data.fDPhiX0(j,0);
        matB(1,2*j+1) = data.fDPhiX0(j,1);
        matB(2,2*j  ) = data.fDPhiX0(j,1);
        matB(2,2*j+1) = data.fDPhiX0(j,0);
    }

    auto force = fForceFunction;
    VecDouble forcingF(fDimension);
    forcingF.setZero();
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);
    
    VecDouble stress(3);

    // TensorToVoigt(Stress,stress);
    stress[0] = Stress.fXX();
    stress[1] = Stress.fYY();
    stress[2] = Stress.fXY();

    // Rhs -= matB.transpose() * fElasticModel->ConstitutiveMatrix() * strain * WJ;
    Rhs -= matB.transpose() * stress * WJ;

    for (int i = nphi; i--; ){
        double shapeFi = data.fPhi[i];
        //External force
        double Fx = forcingF[0] * shapeFi;
        double Fy = forcingF[1] * shapeFi;
        Rhs[2*i  ] += Fx * WJ;
        Rhs[2*i+1] += Fy * WJ;
    };
     
    
};


int PlasticityModel::VariableIndex(const std::string &name) const{
    auto n = fElasticModel->VariableIndex(name);
    if (n != -1){
        return n;
    } else {
        if(!strcmp("PlasticStrain",name.c_str()))    return 101;
    }
    return -1;
    // return fElasticModel->VariableIndex(name);
};

int PlasticityModel::NSolutionVariables(int var) const {
    int val = fElasticModel->NSolutionVariables(var);
    if (val != -1){
        return val;
    } else {
        switch (var)
        {
        case 101:
            return 1;
        
        default:
            return -1;
        }
    }
    return -1;
};

void PlasticityModel::Solution(IntPointData &data, int var, VecDouble &Sol){
    fElasticModel->Solution(data,var,Sol);
    if (var == 101){
        if (data.fPlasticStrain.norm() > 0)
        Sol[0] = data.fPlasticStrain.mean();
        return;
    };
};