#include "ElElasticityPositional2D.h"


void ElElasticityPositional2D::ComputeStiffness(int &index, MatrixDouble &Stiffness){

    double WJ = fIntegData.fWeight * fIntegData.fJacA0 * getIntegPointWeightFunction(index);
    int DIM = Mesh()->Dimension();
    if (DIM != 2) PanicButton();
    
    
    ComputeCurrentJacobian(index);
    ComputeCurrentSpatialDerivatives();

    auto dphi_dx = fIntegData.fDPhiX0;

    //COMPUTE A1
    auto dx_dxsi = fIntegData.fA0;
    auto dy_dxsi = fIntegData.fA1;
    auto j0 = fIntegData.fJacA0;

    //dy_dx
    auto dy_dx = dy_dxsi * fIntegData.fA0Inv.transpose();
    //jacobian
    auto jac = dy_dx.determinant();

    //Green-Lagrange strain tensor
    MatrixDouble E(DIM,DIM);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            E(i,j) = 0.5 * (dy_dx(0,i) * dy_dx(0,j) + dy_dx(1,i) * dy_dx(1,j));
    E(0,0) -= 0.5; E(1,1) -= 0.5;

    //Second Piola-Kirchhoff stress tensor
    MatrixDouble S(DIM,DIM);
    S.setZero();
    double young = Mesh()->getProblemParameters().GetYoungModulus();
    double poisson = Mesh()->getProblemParameters().GetPoissonRatio();
    S(0,0) = young / (1.0-(poisson*poisson)) * (E(0,0) + poisson * E(1,1));
    S(0,1) = 2.0 * (young / (2.0 * (1.0+poisson))) * E(0,1);
    S(1,0) = S(0,1);
    S(1,1) = young / (1.0-(poisson*poisson)) * (E(1,1) + poisson * E(0,0));
    //     S[0][0] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[0][0] + poisson * E[1][1]);
    //     S[1][1] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[1][1] + poisson * E[0][0]);


    //element rhs vector
    for (int a = 0; a < fMesh->NElNodes(); a++){
        for (int k = 0; k < DIM; k++){
            MatrixDouble dE_dyak(DIM,DIM);
            dE_dyak.setZero();
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++)
                    dE_dyak(i,j) = 0.5 * (dphi_dx(a,i) * dy_dx(k,j) + dy_dx(k,i) * dphi_dx(a,j));
            
            //internal force
            // double f = 0.0;
            // for (int i = 0; i < DIM; i++)
            //     for (int j = 0; j < DIM; j++)
            //         f += S(i,j) * dE_dyak(i,j);

            // rhsVector[2 * a + k] -= f * fIntegData.fWeight * j0;

            //element tangent matrix
            for (int b = 0; b < fMesh->NElNodes(); b++){
                for (int l = 0; l < DIM; l++){
                    MatrixDouble dE_dybl(DIM,DIM);
                    dE_dybl.setZero();
                    for (int i = 0; i < DIM; i++)
                        for (int j = 0; j < DIM; j++)
                            dE_dybl(i,j) = 0.5 * (dphi_dx(b,i) * dy_dx(l,j) + dy_dx(l,i) * dphi_dx(b,j));

                    MatrixDouble d2E_dyakbl(DIM,DIM);
                    d2E_dyakbl.setZero();
                    if (k==l)
                        for (int i = 0; i < DIM; i++)
                            for (int j = 0; j < DIM; j++)
                                d2E_dyakbl(i,j) = 0.5 * (dphi_dx(a,i) * dphi_dx(b,j) + dphi_dx(b,i) * dphi_dx(a,j));
                    else
                        for (int i = 0; i < DIM; i++)
                            for (int j = 0; j < DIM; j++)
                                d2E_dyakbl(i,j) = 0.0; 

                    MatrixDouble dS_dybl(DIM,DIM);
                    dS_dybl.setZero();
                    dS_dybl(0,0) = young / (1.0-(poisson*poisson)) * (dE_dybl(0,0) + poisson * dE_dybl(1,1));
                    dS_dybl(0,1) = 2.0 * (young / (2.0 * (1.0+poisson))) * dE_dybl(0,1);
                    dS_dybl(1,0) = 2.0 * (young / (2.0 * (1.0+poisson))) * dE_dybl(1,0);
                    dS_dybl(1,1) = young / (1.0-(poisson*poisson)) * (dE_dybl(1,1) + poisson * dE_dybl(0,0));

                    //elastic and geometric componentes of tangent matrix
                    double e = 0.0;
                    for (int i = 0; i < DIM; i++)
                        for (int j = 0; j < DIM; j++)
                            e += dS_dybl(i,j) * dE_dyak(i,j) + S(i,j) * d2E_dyakbl(i,j);

                    //mass matrix
                    double m = 0.;
                    // (k==l)? m = (1.0 / (0.25 * dTime_)) *0* fMesh->getNumericalIntegration()-> phi_(a,index) * fMesh->getNumericalIntegration()-> phi_(b,index) : m = 0.0;

                    Stiffness(2 * a + k,2 * b + l) += (e+m) * j0 * fIntegData.fWeight;
                }
            }
        }
    }


}

void ElElasticityPositional2D::ComputeResidual(int &index, VecDouble &Rhs){

    auto dphi_dx = fIntegData.fDPhiX0;

    //COMPUTE A1
    auto dx_dxsi = fIntegData.fA0;
    auto dy_dxsi = fIntegData.fA1;
    auto j0 = fIntegData.fJacA0;

    //dy_dx
    auto dy_dx = dy_dxsi * fIntegData.fA0Inv.transpose();
    //jacobian
    auto jac = dy_dx.determinant();

    //Green-Lagrange strain tensor
    MatrixDouble E(DIM,DIM);
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            E(i,j) = 0.5 * (dy_dx(0,i) * dy_dx(0,j) + dy_dx(1,i) * dy_dx(1,j));
    E(0,0) -= 0.5; E(1,1) -= 0.5;

    //Second Piola-Kirchhoff stress tensor
    MatrixDouble S(DIM,DIM);
    S.setZero();
    double young = Mesh()->getProblemParameters().GetYoungModulus();
    double poisson = Mesh()->getProblemParameters().GetPoissonRatio();
    S(0,0) = young / (1.0-(poisson*poisson)) * (E(0,0) + poisson * E(1,1));
    S(0,1) = 2.0 * (young / (2.0 * (1.0+poisson))) * E(0,1);
    S(1,0) = S(0,1);
    S(1,1) = young / (1.0-(poisson*poisson)) * (E(1,1) + poisson * E(0,0));
    //     S[0][0] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[0][0] + poisson * E[1][1]);
    //     S[1][1] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[1][1] + poisson * E[0][0]);


    //element rhs vector
    for (int a = 0; a < fMesh->NElNodes(); a++){
        for (int k = 0; k < DIM; k++){
            MatrixDouble dE_dyak(DIM,DIM);
            dE_dyak.setZero();
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++)
                    dE_dyak(i,j) = 0.5 * (dphi_dx(a,i) * dy_dx(k,j) + dy_dx(k,i) * dphi_dx(a,j));
            
            //internal force
            double f = 0.0;
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++)
                    f += S(i,j) * dE_dyak(i,j);

            Rhs[2 * a + k] -= f * fIntegData.fWeight * j0;
        }
    }
};

void ElElasticityPositional2D::ComputeError(VecDouble &errors){
    std::cout << "Not implemented yet\n";
    return;
    // PanicButton();
}

void ElElasticityPositional2D::ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs){

    for (int i = Mesh()->NElNodes(); i--; ){
        int nstate = Mesh()->NodeVec()[getConnectivity()[i]]->GetNStateVariables();
        for (int istate = 0; istate < nstate; istate++){
            if ((Mesh()->NodeVec()[getConnectivity()[i]] -> getConstrains(istate) == 1) ||
                (Mesh()->NodeVec()[getConnectivity()[i]] -> getConstrains(istate) == 3))  {
                for (int j = Mesh()->NElNodes()*nstate; j--; ){
                    Stiffness(nstate*i+istate,j) = 0.;
                    Stiffness(j,nstate*i+istate) = 0.;
                };
                Stiffness(nstate*i+istate,nstate*i+istate) = 1.;
                Rhs[nstate*i+istate] = 0.;
            }
        }
    }

}