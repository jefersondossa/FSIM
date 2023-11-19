#include "LinearHardening.h"
#include "ElasticTruss.h"
#include "Elasticity2D.h"

LinearHardening::LinearHardening(WeakForm *elast, double hardModulus, double yield) : PlasticityModel(){
    fHardening = hardModulus;
    fYield = yield;
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

    ElasticTruss *truss = dynamic_cast<ElasticTruss* >(fElasticModel);
    Elasticity2D *mat2d = dynamic_cast<Elasticity2D* >(fElasticModel);
    if (truss) {
        fConstitutiveMatrix = truss->ConstitutiveMatrix();
        fYoungModulus = truss->YoungModulus();
    }
    if (mat2d) {
        fConstitutiveMatrix = mat2d->ConstitutiveMatrix();
        fYoungModulus = mat2d->YoungModulus();
    }
}

void LinearHardening::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    fElasticModel->ComputeStiffness(index,data,Stiffness);
};
    
void LinearHardening::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){
    fElasticModel->ComputeResidual(index,data,Rhs);
};
    
void LinearHardening::ComputeError(IntPointData &data, VecDouble &errors){
    fElasticModel->ComputeError(data,errors);

};

int LinearHardening::VariableIndex(const std::string &name) const{
    return fElasticModel->VariableIndex(name);

};

int LinearHardening::NSolutionVariables(int var) const {
    return fElasticModel->NSolutionVariables(var);
};

void LinearHardening::Solution(IntPointData &data, int var, VecDouble &Sol){
    fElasticModel->Solution(data,var,Sol);
};

void LinearHardening::ComputePlasticStrain(IntPointData &data, MatrixDouble &plasticstrain, MatrixDouble &totalstrain){
    //Check if the Integration Point is in the elastic region
    //1 - Compute integration point stress
    int var = fElasticModel->VariableIndex("Stress");
    int nsol = fElasticModel->NSolutionVariables(var);
    VecDouble Sol(nsol);
    MatrixDouble StressTensor(fRealDimension,fRealDimension);
    StressTensor.setZero();
    fElasticModel->Solution(data,var,Sol);
    //Compute the principal stress'
    BuildStressTensor(Sol,StressTensor);
    VecDouble PrincipalStress(fRealDimension);
    ComputePrincipalStress(StressTensor,PrincipalStress);


    if (PrincipalStress[0] > fYield){
        MatrixDouble DeltaStrain = fConstitutiveMatrix.inverse() * StressTensor - totalstrain;
        MatrixDouble DeltaPlasticStrain = (fYoungModulus/(fYoungModulus+fHardening)) * DeltaStrain;
        MatrixDouble DeltaSigma = (fYoungModulus*fHardening/(fYoungModulus+fHardening)) * DeltaStrain;
        std::cout << "Need to do something\n";
    } else {
        totalstrain = fConstitutiveMatrix.inverse() * StressTensor; 
    }

}

