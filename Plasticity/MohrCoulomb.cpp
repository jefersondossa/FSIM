#include "MohrCoulomb.h"
#include "ElasticTruss.h"
#include "Elasticity2D.h"

MohrCoulomb::MohrCoulomb(WeakForm *elast, double cohesion, double intfriction) : PlasticityModel(elast){
    fCohesion = cohesion;
    fInternalFriction = intfriction;

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

void MohrCoulomb::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    fElasticModel->ComputeStiffness(index,data,Stiffness);
    // Change stiffness;
    if (fPlasticStrain.rows()>0){
        PanicButton();
        // if (fRealDimension == 1){
        //     Stiffness *= fHardening/(fYoungModulus+fHardening);
        // } else {
        //     PanicButton();
        // }
    }
};
    
void MohrCoulomb::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){
    fElasticModel->ComputeResidual(index,data,Rhs);
    if (fPlasticStrain.rows()>0){
        PanicButton();
        // if (fRealDimension == 1){
        //     int nphi = data.fPhi.size();
        //     double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
        //     double elementLenght = 2.*data.fJacA0;
        //     double K = fYoungModulus / elementLenght;

        //     MatrixDouble rotation(fDimension*nphi,fDimension*nphi);
        //     MatrixDouble matB(fDimension,fDimension*nphi);
        //     rotation.setZero();
        //     matB.setZero();
        //     double cosa = data.fAxes(0,0) / data.fJacA0;
        //     double sina = data.fAxes(1,0) / data.fJacA0;
        //     double check = sina*sina+cosa*cosa;
        //     for (int j = 0; j < nphi; j++){
        //         // for (int i = 0; i < fDimension; i++){
        //             matB(0,fDimension*j) = data.fDPhiX0(j,0);
        //         // }
        //         rotation(2*j  ,2*j  ) = cosa;
        //         rotation(2*j+1,2*j  ) = sina;
        //         rotation(2*j  ,2*j+1) = -sina;
        //         rotation(2*j+1,2*j+1) = cosa;
        //     }

        //     VecDouble sol(2);// = data.fSol;
        //     sol[0] = -fPlasticStrain(0,0)*sina;
        //     sol[1] = +fPlasticStrain(0,0)*cosa;

        //     // std::cout << "rotation =\n"<< rotation << std::endl;
        //     Rhs += rotation * matB.transpose() * sol * WJ * elementLenght * K;
        // } else {
        //     PanicButton();
        // }
    }
    
};
    
void MohrCoulomb::ComputeError(IntPointData &data, VecDouble &errors){
    fElasticModel->ComputeError(data,errors);

};

int MohrCoulomb::VariableIndex(const std::string &name) const{
    return fElasticModel->VariableIndex(name);

};

int MohrCoulomb::NSolutionVariables(int var) const {
    return fElasticModel->NSolutionVariables(var);
};

void MohrCoulomb::Solution(IntPointData &data, int var, VecDouble &Sol){
    fElasticModel->Solution(data,var,Sol);
};

void MohrCoulomb::ComputePlasticStrain(IntPointData &data, MatrixDouble &plasticstrain, MatrixDouble &totalstrain){
    //Check if the Integration Point is in the elastic region
    //1 - Compute integration point stress
    int var = fElasticModel->VariableIndex("Stress");
    int nsol = fElasticModel->NSolutionVariables(var);
    VecDouble Sol(nsol);
    MatrixDouble StressTensor(fRealDimension,fRealDimension);
    StressTensor.setZero();
    fElasticModel->Solution(data,var,Sol);
    //Compute the principal stress'
    VoigtToTensor(StressTensor,Sol);
    VecDouble PrincipalStress(fRealDimension);
    ComputePrincipalStress(StressTensor,PrincipalStress);

    double criterion = fYoungModulus;
    if (fabs(PrincipalStress[0]) > criterion){
        std::cout << "Please Implement me\n";
    } else {
        totalstrain = fConstitutiveMatrix.inverse() * StressTensor;
        fPlasticStrain.resize(0,0);
        fTotalStrain.resize(0,0);
    }

}