#include "ElasticityPositional2D.h"

ElasticityPositional2D::ElasticityPositional2D(int matid, double young, double poisson, bool planes) : WeakForm() {
    this->fMatId = matid;
    fDimension = 2;
    fNState = 2;
    fYoungModulus = young;
    fPoissonRatio = poisson;
    fPlaneStress = planes;
};

void ElasticityPositional2D::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    if (!data.fNeedsDSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
    }
    int nphi = data.fPhi.size();
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    
    auto dphi_dx = data.fDPhiX0;

    //COMPUTE A1
    auto dx_dxsi = data.fA0;
    auto dy_dxsi = data.fA1;
    auto j0 = data.fJacA0;

    //dy_dx
    auto dy_dx = dy_dxsi * data.fA0Inv;
    //jacobian
    auto jac = dy_dx.determinant();

    //Green-Lagrange strain tensor
    MatrixDouble E(fDimension,fDimension);
    for (int i = 0; i < fDimension; i++)
        for (int j = 0; j < fDimension; j++)
            E(i,j) = 0.5 * (dy_dx(0,i) * dy_dx(0,j) + dy_dx(1,i) * dy_dx(1,j));
    E(0,0) -= 0.5; E(1,1) -= 0.5;

    //Second Piola-Kirchhoff stress tensor
    MatrixDouble S(fDimension,fDimension);
    S.setZero();
    if (fPlaneStress){
        S(0,0) = fYoungModulus / (1.0-(fPoissonRatio*fPoissonRatio)) * (E(0,0) + fPoissonRatio * E(1,1));
        S(0,1) = 2.0 * (fYoungModulus / (2.0 * (1.0+fPoissonRatio))) * E(0,1);
        S(1,0) = S(0,1);
        S(1,1) = fYoungModulus / (1.0-(fPoissonRatio*fPoissonRatio)) * (E(1,1) + fPoissonRatio * E(0,0));
    } else {
        std::cout << "Pleas implement me " << std::endl;
        PanicButton();
    }
    //     S[0][0] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[0][0] + poisson * E[1][1]);
    //     S[1][1] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[1][1] + poisson * E[0][0]);


    //element rhs vector
    for (int a = 0; a < nphi; a++){
        for (int k = 0; k < fDimension; k++){
            MatrixDouble dE_dyak(fDimension,fDimension);
            dE_dyak.setZero();
            for (int i = 0; i < fDimension; i++)
                for (int j = 0; j < fDimension; j++)
                    dE_dyak(i,j) = 0.5 * (dphi_dx(a,i) * dy_dx(k,j) + dy_dx(k,i) * dphi_dx(a,j));
            
            //internal force
            // double f = 0.0;
            // for (int i = 0; i < fDimension; i++)
            //     for (int j = 0; j < fDimension; j++)
            //         f += S(i,j) * dE_dyak(i,j);

            // rhsVector[2 * a + k] -= f * data.fWeight * j0;

            //element tangent matrix
            for (int b = 0; b < nphi; b++){
                for (int l = 0; l < fDimension; l++){
                    MatrixDouble dE_dybl(fDimension,fDimension);
                    dE_dybl.setZero();
                    for (int i = 0; i < fDimension; i++)
                        for (int j = 0; j < fDimension; j++)
                            dE_dybl(i,j) = 0.5 * (dphi_dx(b,i) * dy_dx(l,j) + dy_dx(l,i) * dphi_dx(b,j));

                    MatrixDouble d2E_dyakbl(fDimension,fDimension);
                    d2E_dyakbl.setZero();
                    if (k==l)
                        for (int i = 0; i < fDimension; i++)
                            for (int j = 0; j < fDimension; j++)
                                d2E_dyakbl(i,j) = 0.5 * (dphi_dx(a,i) * dphi_dx(b,j) + dphi_dx(b,i) * dphi_dx(a,j));
                    else
                        for (int i = 0; i < fDimension; i++)
                            for (int j = 0; j < fDimension; j++)
                                d2E_dyakbl(i,j) = 0.0; 

                    MatrixDouble dS_dybl(fDimension,fDimension);
                    dS_dybl.setZero();
                    dS_dybl(0,0) = fYoungModulus / (1.0-(fPoissonRatio*fPoissonRatio)) * (dE_dybl(0,0) + fPoissonRatio * dE_dybl(1,1));
                    dS_dybl(0,1) = 2.0 * (fYoungModulus / (2.0 * (1.0+fPoissonRatio))) * dE_dybl(0,1);
                    dS_dybl(1,0) = 2.0 * (fYoungModulus / (2.0 * (1.0+fPoissonRatio))) * dE_dybl(1,0);
                    dS_dybl(1,1) = fYoungModulus / (1.0-(fPoissonRatio*fPoissonRatio)) * (dE_dybl(1,1) + fPoissonRatio * dE_dybl(0,0));

                    //elastic and geometric componentes of tangent matrix
                    double e = 0.0;
                    for (int i = 0; i < fDimension; i++)
                        for (int j = 0; j < fDimension; j++)
                            e += dS_dybl(i,j) * dE_dyak(i,j) + S(i,j) * d2E_dyakbl(i,j);

                    //mass matrix
                    double m = 0.;
                    // (k==l)? m = (1.0 / (0.25 * dTime_)) *0* data.fPhi[a] * data.fPhi[b] : m = 0.0;

                    Stiffness(2 * a + k,2 * b + l) += (e+m) * j0 * data.fWeight;
                }
            }
        }
    }


}

void ElasticityPositional2D::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    auto dphi_dx = data.fDPhiX0;
    int nphi = data.fPhi.size();
    //COMPUTE A1
    auto dx_dxsi = data.fA0;
    auto dy_dxsi = data.fA1;
    auto j0 = data.fJacA0;

    //dy_dx
    auto dy_dx = dy_dxsi * data.fA0Inv;
    //jacobian
    auto jac = data.fJacA1;

    //Green-Lagrange strain tensor
    MatrixDouble E(fDimension,fDimension);
    for (int i = 0; i < fDimension; i++)
        for (int j = 0; j < fDimension; j++)
            E(i,j) = 0.5 * (dy_dx(0,i) * dy_dx(0,j) + dy_dx(1,i) * dy_dx(1,j));
    E(0,0) -= 0.5; E(1,1) -= 0.5;

    //Second Piola-Kirchhoff stress tensor
    MatrixDouble S(fDimension,fDimension);
    S.setZero();
    double young = fYoungModulus;
    double poisson = fPoissonRatio;
    if (fPlaneStress){
        S(0,0) = young / (1.0-(poisson*poisson)) * (E(0,0) + poisson * E(1,1));
        S(0,1) = 2.0 * (young / (2.0 * (1.0+poisson))) * E(0,1);
        S(1,0) = S(0,1);
        S(1,1) = young / (1.0-(poisson*poisson)) * (E(1,1) + poisson * E(0,0));
    } else {
        std::cout << "Pleas implement me " << std::endl;
        PanicButton();
    }
    
    //     S[0][0] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[0][0] + poisson * E[1][1]);
    //     S[1][1] = (young / ( (1.0+poisson) * (1.0-2.0*poisson))) * ((1.0-poisson) * E[1][1] + poisson * E[0][0]);


    //element rhs vector
    for (int a = 0; a < nphi; a++){
        for (int k = 0; k < fDimension; k++){
            MatrixDouble dE_dyak(fDimension,fDimension);
            dE_dyak.setZero();
            for (int i = 0; i < fDimension; i++)
                for (int j = 0; j < fDimension; j++)
                    dE_dyak(i,j) = 0.5 * (dphi_dx(a,i) * dy_dx(k,j) + dy_dx(k,i) * dphi_dx(a,j));
            
            //internal force
            double f = 0.0;
            for (int i = 0; i < fDimension; i++)
                for (int j = 0; j < fDimension; j++)
                    f += S(i,j) * dE_dyak(i,j);

            Rhs[2 * a + k] -= f * data.fWeight * j0;
        }
    }
};

void ElasticityPositional2D::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();

}