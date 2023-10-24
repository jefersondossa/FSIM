#include "ElElasticityPositional2D.h"


void ElElasticityPositional2D::ComputeStiffness(int &index, MatrixDouble &Stiffness){

    double WJ = fIntegData.fWeight * fIntegData.fJacA0 * getIntegPointWeightFunction(index);
    int DIM = Mesh()->Dimension();
    if (DIM != 2) PanicButton();

    //COMPUTE A1
    MatrixDouble dy_dxsi(DIM,DIM), dx_dxsi(DIM,DIM); //row = cartesian, column = parametric
    dy_dxsi.setZero();
    dx_dxsi.setZero();
    
    for (int i = 0; i < Mesh()->NElNodes(); i++){
        VecDouble currentCoord(2);
        currentCoord[0] = Mesh()->NodeVec()[this->fConnect[i]] -> GetSolution(0);
        currentCoord[1] = Mesh()->NodeVec()[this->fConnect[i]] -> GetSolution(1);

        for (int k = 0; k < DIM; k++)
            for (int l = 0; l < DIM; l++)
                dy_dxsi(k,l) += currentCoord[k] * fMesh->getNumericalIntegration()->dphi_[i](l,index);
    }

    //dy_dx
    // double dy_dx(DIM,DIM,0.);
    // for (int i = 0; i < DIM; i++)
    //     for (int j = 0; j < DIM; j++)
    //         dy_dx(i,j) = dy_dxsi(i,0) * ainv(0,j) + dy_dxsi(i,1) * dxsi_dx(1,j);

    // //jacobian
    // double jac = dy_dx[0][0] * dy_dx[1][1] - dy_dx[0][1] * dy_dx[1][0];

    // //Green-Lagrange strain tensor
    // double E[DIM][DIM] = {};
    // for (int i = 0; i < DIM; i++)
    //     for (int j = 0; j < DIM; j++)
    //         E[i][j] = 0.5 * (dy_dx[0][i] * dy_dx[0][j] + dy_dx[1][i] * dy_dx[1][j]);
    // E[0][0] -= 0.5; E[1][1] -= 0.5;

    // //Second Piola-Kirchhoff stress tensor
    // double S[DIM][DIM] = {};
    // double young = 1. / j0;
    // double poisson = 0.3;
    // S[0][0] = young / (1.0-(poisson*poisson)) * (E[0][0] + poisson * E[1][1]);
    // S[0][1] = 2.0 * (young / (2.0 * (1.0+poisson))) * E[0][1];
    // S[1][0] = S[0][1];
    // S[1][1] = young / (1.0-(poisson*poisson)) * (E[1][1] + poisson * E[0][0]);
    // //     S[0][0] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[0][0] + poisson * E[1][1]);
    // //     S[1][1] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[1][1] + poisson * E[0][0]);


    // //element rhs vector
    // for (int a = 0; a < fMesh->NElNodes(); a++){
    //     for (int k = 0; k < DIM; k++){
    //         double dE_dyak[DIM][DIM];
    //         for (int i = 0; i < DIM; i++)
    //             for (int j = 0; j < DIM; j++)
    //                 dE_dyak[i][j] = 0.5 * (dphi_dx(a,i) * dy_dx[k][j] + dy_dx[k][i] * dphi_dx(a,j));
            
    //         //internal force
    //         double f = 0.0;
    //         for (int i = 0; i < DIM; i++)
    //             for (int j = 0; j < DIM; j++)
    //                 f += S[i][j] * dE_dyak[i][j];

    //         double vel = 0.0;
    //         for (int i = 0; i < fMesh->NElNodes(); i++)
    //             vel += fMesh->getNumericalIntegration()-> phi_(i,index) * fMesh->NodeVec()[fConnect[i]]->getMeshVelocity(k);

    //         double c =  fMesh->getNumericalIntegration()-> phi_(a,index) * vel*0;

    //         //domain force
    //         // double b;
    //         // (k==1) ? b = phi(a) * density * gravity : 0.0;

    //         rhsVector[2 * a + k] -= (f+c) * weight_ * j0;

    //         //element tangent matrix
    //         for (int b = 0; b < fMesh->NElNodes(); b++){
    //             for (int l = 0; l < DIM; l++){
    //                 double dE_dybl[DIM][DIM];
    //                 for (int i = 0; i < DIM; i++)
    //                     for (int j = 0; j < DIM; j++)
    //                         dE_dybl[i][j] = 0.5 * (dphi_dx(b,i) * dy_dx[l][j] + dy_dx[l][i] * dphi_dx(b,j));

    //                 double d2E_dyakbl[DIM][DIM];
    //                 if (k==l)
    //                     for (int i = 0; i < DIM; i++)
    //                         for (int j = 0; j < DIM; j++)
    //                             d2E_dyakbl[i][j] = 0.5 * (dphi_dx(a,i) * dphi_dx(b,j) + dphi_dx(b,i) * dphi_dx(a,j));
    //                 else
    //                     for (int i = 0; i < DIM; i++)
    //                         for (int j = 0; j < DIM; j++)
    //                             d2E_dyakbl[i][j] = 0.0; 

    //                 double dS_dybl[DIM][DIM];
    //                 dS_dybl[0][0] = young / (1.0-(poisson*poisson)) * (dE_dybl[0][0] + poisson * dE_dybl[1][1]);
    //                 dS_dybl[0][1] = 2.0 * (young / (2.0 * (1.0+poisson))) * dE_dybl[0][1];
    //                 dS_dybl[1][0] = 2.0 * (young / (2.0 * (1.0+poisson))) * dE_dybl[1][0];
    //                 dS_dybl[1][1] = young / (1.0-(poisson*poisson)) * (dE_dybl[1][1] + poisson * dE_dybl[0][0]);

    //                 //elastic and geometric componentes of tangent matrix
    //                 double e = 0.0;
    //                 for (int i = 0; i < DIM; i++)
    //                     for (int j = 0; j < DIM; j++)
    //                         e += dS_dybl[i][j] * dE_dyak[i][j] + S[i][j] * d2E_dyakbl[i][j];

    //                 //mass matrix
    //                 double m;
    //                 (k==l)? m = (1.0 / (0.25 * dTime_)) *0* fMesh->getNumericalIntegration()-> phi_(a,index) * fMesh->getNumericalIntegration()-> phi_(b,index) : m = 0.0;

    //                 jacobianNRMatrix(2 * a + k,2 * b + l) += (e+m) * j0 * weight_;
    //             }
    //         }
    //     }
    // }


}

void ElElasticityPositional2D::ComputeResidual(int &index, VecDouble &Rhs){

    VecDouble fieldForce = Mesh()->getProblemParameters().GetFieldForce();
    auto force = Mesh()->getProblemParameters().getForcingFunction();
    int dim = Mesh()->Dimension();

    //Velocity Derivatives
    MatrixDouble du_dx(dim,dim), duprev_dx(dim,dim), duna_dx(dim,dim);
    interpolateSolDerivatives(du_dx);
    duna_dx = du_dx;

    double WJ = fIntegData.fWeight * fIntegData.fJacA0  * getIntegPointWeightFunction(index);

    VecDouble forcingF(1);
    VecDouble x_ = getIntegPointCoordinatesValue(index);
    if (force) force(x_,forcingF);

    for (int i = Mesh()->NElNodes(); i--; ){
        double shapeFi = Mesh()->getNumericalIntegration()-> phi_(i,index);

        //Viscosity
        double K = 0.;
        for (int l=dim; l--; ) K += fIntegData.fDPhiX0(i,l) * duna_dx(0,l);

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