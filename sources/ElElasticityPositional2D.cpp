#include "ElElasticityPositional2D.h"


void ElElasticityPositional2D::ComputeStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, MatrixDouble &Stiffness){

    double WJ = weight_ * djac_ * getIntegPointWeightFunction(index);

    for (int i = Mesh()->NElNodes(); i-- ; ){       
        for (int j = Mesh()->NElNodes(); j-- ; ){            
            for (int k = Mesh()->Dimension(); k--;  ){
                //Diffusion matrix
                double K = dphi_dx(i,k) * dphi_dx(j,k);
                Stiffness(i,j) += K * WJ;
            }
        };
    };


}

void ElElasticityPositional2D::ComputeResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, VecDouble &Rhs){

    VecDouble fieldForce = Mesh()->getProblemParameters().GetFieldForce();
    auto force = Mesh()->getProblemParameters().getForcingFunction();
    int dim = Mesh()->Dimension();

    //Velocity Derivatives
    MatrixDouble du_dx(dim,dim), duprev_dx(dim,dim), duna_dx(dim,dim);
    interpolateVelDerivatives(dphi_dx, du_dx, duprev_dx);
    duna_dx = du_dx;

    double WJ = weight_ * djac_  * getIntegPointWeightFunction(index);

    VecDouble forcingF(1);
    VecDouble x_ = getIntegPointCoordinatesValue(index);
    if (force) force(x_,forcingF);

    for (int i = Mesh()->NElNodes(); i--; ){
        double shapeFi = Mesh()->getNumericalIntegration()-> phi_(i,index);

        //Viscosity
        double K = 0.;
        for (int l=dim; l--; ) K += dphi_dx(i,l) * duna_dx(0,l);

        //External force
        double F = (fieldForce[0] + forcingF[0]) * shapeFi;
        Rhs[i] += (-K + F) * WJ;
    };
};

void ElElasticityPositional2D::ComputeError(VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();
}

void ElElasticityPositional2D::ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs){

    for (int i = Mesh()->NElNodes(); i--; ){
        if ((Mesh()->NodeVec()[getConnectivity()[i]] -> getConstrains(0) == 1) ||
            (Mesh()->NodeVec()[getConnectivity()[i]] -> getConstrains(0) == 3))  {
            for (int j = Mesh()->NElNodes(); j--; ){
                Stiffness(i,j) = 0.;
                Stiffness(j,i) = 0.;
            };
            Stiffness(i,i) = 1.;
            Rhs[i] = 0.;
        }
    }

}