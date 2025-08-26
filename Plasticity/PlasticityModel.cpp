
#include "PlasticityModel.h"
#include "ElasticTruss.h"
#include "Elasticity2D.h"
#include "Elasticity3D.h"
#include "PositionalTruss.h"
#include "ElasticityPositional2D.h"

PlasticityModel::PlasticityModel(WeakForm *elast){
    fElasticModel = elast;
    fNState = elast->NState();
    fMatId = elast->Id();
    int var = fElasticModel->VariableIndex("Stress");
    int nsol = fElasticModel->NSolutionVariables(var);
    fRealDimension = elast->Dimension();
    // if (nsol == 1) fRealDimension = 1;
    // if (nsol == 3) fRealDimension = 2;
    // if (nsol == 6) fRealDimension = 3;
    fNStressComponents = nsol;
    fDimension = elast->Dimension();

    ElasticTruss *truss = dynamic_cast<ElasticTruss* >(fElasticModel);
    Elasticity2D *mat2d = dynamic_cast<Elasticity2D* >(fElasticModel);
    Elasticity3D *mat3d = dynamic_cast<Elasticity3D* >(fElasticModel);
    PositionalTruss *postruss = dynamic_cast<PositionalTruss* >(fElasticModel);
    ElasticityPositional2D *pos2d = dynamic_cast<ElasticityPositional2D* >(fElasticModel);
    if (truss) {
        fYoungModulus = truss->YoungModulus();
        fPoissonRatio = 0.;
    }
    if (mat2d) {
        fYoungModulus = mat2d->YoungModulus();
        fPoissonRatio = mat2d->PoissonRatio();
        fPlaneStress = mat2d->PlaneState();
        fThickness = mat2d->Thickness();
    }
    if (mat3d) {
        fYoungModulus = mat3d->YoungModulus();
        fPoissonRatio = mat3d->PoissonRatio();
        fPlaneStress = false;
    }
    if (postruss) {
        fYoungModulus = postruss->YoungModulus();
        fPoissonRatio = 0.;
    }
    if (pos2d) {
        fYoungModulus = pos2d->YoungModulus();
        fPoissonRatio = pos2d->PoissonRatio();
        fPlaneStress = pos2d->PlaneState();
    }

    Tensor3D Id2;
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
    
    //Transform 4th order Tensor3D into matrix: https://wiki.seg.org/wiki/Voigt_notation
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

void PlasticityModel::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs, Tensor3D &Stress){

    // fElasticModel->ComputeResidual(index,data,Rhs);
    // return;
    int nphi = data.fPhi.size();

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index] * fThickness;
    MatrixDouble matB;

    if (this->Dimension() == 2){
        matB.resize(3,2*nphi);
        matB.setZero();
        for (int j = 0; j < nphi; j++){
            matB(0,2*j  ) = data.fDPhiX0(0,j);
            matB(1,2*j+1) = data.fDPhiX0(1,j);
            matB(2,2*j  ) = data.fDPhiX0(1,j);
            matB(2,2*j+1) = data.fDPhiX0(0,j);
        }
    } else if (this->Dimension() == 3){
        matB.resize(6,3*nphi);
        matB.setZero();
        for (int j = 0; j < nphi; j++){
            matB(0,3*j  ) = data.fDPhiX0(0,j);
            matB(1,3*j+1) = data.fDPhiX0(1,j);
            matB(2,3*j+2) = data.fDPhiX0(2,j);
            matB(5,3*j  ) = data.fDPhiX0(1,j);
            matB(5,3*j+1) = data.fDPhiX0(0,j);
            matB(3,3*j+1) = data.fDPhiX0(2,j);
            matB(3,3*j+2) = data.fDPhiX0(1,j);
            matB(4,3*j  ) = data.fDPhiX0(2,j);
            matB(4,3*j+2) = data.fDPhiX0(0,j);
        }
    }

    auto force = fForceFunction;
    VecDouble forcingF(fDimension);
    forcingF.setZero();
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);
    
    VecDouble stress;

    if (fDimension == 2){
        stress.resize(3);
        stress[0] = Stress.fXX();
        stress[1] = Stress.fYY();
        stress[2] = Stress.fXY();
    } else if (fDimension == 3){
        stress.resize(6);
        stress[0] = Stress.fXX();
        stress[1] = Stress.fYY();
        stress[2] = Stress.fZZ();
        stress[3] = Stress.fYZ();
        stress[4] = Stress.fXZ();
        stress[5] = Stress.fXY();
    }
    

    // Rhs -= matB.transpose() * fElasticModel->ConstitutiveMatrix() * strain * WJ;
    Rhs -= matB.transpose() * stress * WJ;

    for (int i = nphi; i--; ){
        double shapeFi = data.fPhi[i];
        //External force
        for (int k = 0; k < fDimension; k++){
            Rhs[fDimension*i+k] += forcingF[k] * shapeFi * WJ;
        }
    };
     
    
};


int PlasticityModel::VariableIndex(const std::string &name) const{
    auto n = fElasticModel->VariableIndex(name);
    if (n != -1){
        return n;
    } else {
        if(!strcmp("PlasticStrain",name.c_str()))    return 101;
        if(!strcmp("RealStress",name.c_str()))    return 102;
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
        case 102:
            return 3;
        
        default:
            return -1;
        }
    }
    return -1;
};

void PlasticityModel::Solution(IntPointData &data, int var, VecDouble &Sol){
    fElasticModel->Solution(data,var,Sol);
    if (var == 101){
        if (data.fEffectivePlasticStrain.norm() > 0)
        Sol[0] = data.fEffectivePlasticStrain.mean();
        return;
    };
    if (var == 102){
        MatrixDouble ElasticConstitutive = 2.*fShearModulus*fIdentity4Dev + fBulkModulus*fId2xId2;
        auto stress = data.fElasticStrain[data.fIndex].Multiply(ElasticConstitutive);
        // if (data.fEffectivePlasticStrain.norm() > 0)
        Sol[0] = stress.fXX();
        Sol[1] = stress.fYY();
        Sol[2] = stress.fXY();
        return;
    };
};