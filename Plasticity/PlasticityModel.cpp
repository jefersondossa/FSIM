
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

    fIdentity2.resize(fNStressComponents);
    fIdentity2.setZero();
    fIdentity4.resize(fNStressComponents,fNStressComponents);
    fIdentity4Dev.resize(fNStressComponents,fNStressComponents);
    fIdentity4.setIdentity();
    fIdentity4Dev.setZero();

    if (fRealDimension == 2){
        //Definition of 2nd and 4th order identity tensors  
        fIdentity2[0] = 1.;
        fIdentity2[1] = 1.;
        fIdentity4(2,2) = 0.5;
        fIdentity4Dev = fIdentity4;
        for (int i = 0; i < fRealDimension; i++){
            for (int j = 0; j < fRealDimension; j++){
                fIdentity4Dev(i,j) -= fIdentity2[i] * fIdentity2[j] / 3.; 
            }
        }
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