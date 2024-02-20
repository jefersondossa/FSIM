#include "DruckerPrager.h"

DruckerPrager::DruckerPrager(WeakForm *elast) : PlasticityModel(elast){


}

void DruckerPrager::ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor &Stress){
    
    MatrixDouble fTangentTensor(6,6);
    fTangentTensor.setZero();
    
    std::cout << "Put the tangent elastoplastic tensor here.\n";
    PanicButton();   

    if (fElasticModel->Dimension() == 2){
        VecInt order(3);
        order[0] = XX;
        order[1] = YY;
        order[2] = XY;
        MatrixDouble fTangent2D(3,3);
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                fTangent2D(i,j) = fTangentTensor(order[i],order[j]);
            }
        }
        fElasticModel->ConstitutiveMatrix() = fTangent2D;
    } else if (fElasticModel->Dimension() == 3) {
        fElasticModel->ConstitutiveMatrix() = fTangentTensor;
    } else {
        PanicButton();
    }
    
    fElasticModel->ComputeStiffness(index,data,Stiffness);
};
    
void DruckerPrager::ComputeError(IntPointData &data, VecDouble &errors){
    fElasticModel->ComputeError(data,errors);
};


double DruckerPrager::YieldFunction(int &index, IntPointData &data, Tensor &Stress){
    double YF = 0.;
    std::cout << "Put the yield function here! \n";
    PanicButton();
    if (fPlaneStress){
         
    } else {
        
    }
    return YF;
}

double DruckerPrager::PlasticMultiplier(int &index, IntPointData &data, Tensor &Stress){
    //Newton-Raphson to find plastic multiplier
    std::cout << "Put the plastic multiplier here!\n";
    PanicButton();
    double dGamma = 0.;
    if (fPlaneStress){
        
    } else {
        
    }

    return dGamma;
}

void DruckerPrager::UpdateStateVariables(int &index, IntPointData &data, Tensor &Stress){
    
    std::cout << "Update the state variables here! \n";
    PanicButton();

}   