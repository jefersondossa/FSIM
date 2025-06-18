#include "PositionalFrame2D.h"
#include "IntRule1d.h"


PositionalFrame2D::PositionalFrame2D(int matid, double young, double base, double height) : WeakForm() {
    this->fMatId = matid;
    fDimension = 1;
    fNState = 3;
    fYoungModulus = young;
    fDepth = base;
    fHeight = height;
};


void PositionalFrame2D::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    if (!data.fNeedsDSol || !data.fNeedsSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
        data.fNeedsSol = true;
        data.fSol.resize(fNState);
    }

    IntRule1d ruleEta(4,3);

    //Tangent vector
    VecDouble normalVersor(2), normalVersorUp(2);
    double normTangent = sqrt(data.fAxes0(0,0)*data.fAxes0(0,0) + data.fAxes0(1,0)*data.fAxes0(1,0));
    double normTangentUp = sqrt(data.fAxes1(0,0)*data.fAxes1(0,0) + data.fAxes1(1,0)*data.fAxes1(1,0));
    normalVersor[0] = -data.fAxes0(1,0)/normTangent;
    normalVersor[1] =  data.fAxes0(0,0)/normTangent;
    normalVersorUp[0] = -data.fAxes1(1,0)/normTangentUp;
    normalVersorUp[1] =  data.fAxes1(0,0)/normTangentUp;
    double theta0 = atan2(normalVersor[1],normalVersor[0]);
    // double theta1 = atan2(normalVersorUp[1],normalVersorUp[0]);
    double theta1 = theta0 + data.fSol[2];

    int nphi = data.fPhi.size();

    //Matrix A0 (initial configuration) and A1 (current configuration)
    MatrixDouble A0(2,2), A1(2,2);
    A0.setZero();
    A1.setZero();

    for (int eta = 0; eta < ruleEta.NPoints(); eta++){
        double coordEta = ruleEta.PointList(eta,0);
        double weightEta = ruleEta.WeightList(eta);
        
        ComputeA0A1(data, coordEta, theta0, theta1, A0, A1);
        
        MatrixDouble A0inv = A0.inverse();
        MatrixDouble E;
        ComputeGreenStrain(data, A0inv, A1, E);

        //Jacobian value (Equation 6.50)
        double J0 = A0(0,0) * A0(1,1) - A0(0,1) * A0(1,0);

        //Saint-Venant_Kirchhoff Stress
        //Equation 6.62 to 6.65)
        MatrixDouble S = fYoungModulus * E;
    
        //Terms from the second derivate of matrix A¹ (d²A¹/dtheta_beta dtheta_z)
        for(int alpha = 0; alpha < 3; alpha++){
            for(int beta = 0; beta < nphi; beta++){
                MatrixDouble DA1DY_beta(2,2);
                //Equation 6.69, 6.71 and 6.72
                ComputeDA1_DY(data, alpha, theta1, coordEta, beta, A0inv, A1, DA1DY_beta);
                //6.70 or 6.67
                MatrixDouble DEDy_beta;
                ComputeDE_DY(A0inv, A1, DA1DY_beta, DEDy_beta);

                for (int gama = 0; gama < 3; gama++){
                    for (int zeta = 0; zeta < nphi; zeta++){

                        MatrixDouble DA1DY_zeta(2,2);
                        //Equation 6.69, 6.71 and 6.72
                        ComputeDA1_DY(data, gama, theta1, coordEta, zeta, A0inv, A1, DA1DY_zeta);
                        //6.70 or 6.67
                        MatrixDouble DEDy_zeta;
                        ComputeDE_DY(A0inv, A1, DA1DY_zeta, DEDy_zeta);
                        
                        MatrixDouble D2DA1_DY2(2,2);
                        D2DA1_DY2.setZero();
                        if (alpha==2 && gama==2){
                            //Equation 6.83 to 6.86
                            D2DA1_DY2(0,0) = 0.5 * fHeight * coordEta * (sin(theta1) * data.fPhi[zeta] * data.fPhi[beta] * data.fDSolDAdim(2,0) 
                                                    - cos(theta1) * (data.fPhi[zeta] * data.fDPhi(0,beta) + data.fPhi[beta] * data.fDPhi(0,zeta)));
                            D2DA1_DY2(1,0) =-0.5 * fHeight * coordEta * (cos(theta1) * data.fPhi[zeta] * data.fPhi[beta] * data.fDSolDAdim(2,0) 
                                                    + sin(theta1) * (data.fPhi[zeta] * data.fDPhi(0,beta) + data.fPhi[beta] * data.fDPhi(0,zeta)));
                            D2DA1_DY2(0,1) = -0.5 * fHeight * (cos(theta1) * data.fPhi[zeta] * data.fPhi[beta]); 
                            D2DA1_DY2(1,1) = -0.5 * fHeight * (sin(theta1) * data.fPhi[zeta] * data.fPhi[beta]);

                        }
                        
                        //Equation 6.81
                        MatrixDouble temp1 = A0inv.transpose() * DA1DY_beta.transpose() * DA1DY_zeta* A0inv;
                        MatrixDouble temp2 = A0inv.transpose() * D2DA1_DY2.transpose() * A1 * A0inv;
                        MatrixDouble D2E_DY2 = 0.5 * (temp1 + temp1.transpose() + temp2 + temp2.transpose());

                        //Equation 6.79
                        MatrixDouble DSDy = fYoungModulus * DEDy_zeta;

                        //Equation 6.76
                        double DSDYcDEDy = DoubleContraction(DSDy,DEDy_beta);
                        double ScD2EDY2 = DoubleContraction(S,D2E_DY2);                  

                        int j = 3 * (beta) + alpha;
                        int i = 3 * (zeta) + gama;

                        Stiffness(i,j) += (DSDYcDEDy+ScD2EDY2) * fDepth * weightEta * J0 * data.fWeight;                        
                    }
                }
            }
        } 
    }
}

void PositionalFrame2D::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    auto force = fForceFunction;
    VecDouble forcingF(3);
    forcingF.setZero();
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    IntRule1d ruleEta(4,3);

    //Tangent vector
    VecDouble normalVersor(2);
    double normTangent = sqrt(data.fAxes0(0,0)*data.fAxes0(0,0) + data.fAxes0(1,0)*data.fAxes0(1,0));
    double normTangentUp = sqrt(data.fAxes1(0,0)*data.fAxes1(0,0) + data.fAxes1(1,0)*data.fAxes1(1,0));
    normalVersor[0] = -data.fAxes0(1,0)/normTangent;
    normalVersor[1] =  data.fAxes0(0,0)/normTangent;
    double theta0 = atan2(normalVersor[1],normalVersor[0]);
    // double theta1 = atan2(normalVersorUp[1],normalVersorUp[0]);
    double theta1 = theta0 + data.fSol[2];
    
    int nphi = data.fPhi.size();

    //Matrix A0 (initial configuration) and A1 (current configuration)
    MatrixDouble A0(2,2), A1(2,2);
    A0.setZero();
    A1.setZero();
    
    for (int eta = 0; eta < ruleEta.NPoints(); eta++){
        double coordEta = ruleEta.PointList(eta,0);
        double weightEta = ruleEta.WeightList(eta);
        
        ComputeA0A1(data, coordEta, theta0, theta1, A0, A1);
    
        MatrixDouble A0inv = A0.inverse();

        MatrixDouble E(2,2);
        ComputeGreenStrain(data, A0inv, A1, E);

        //Jacobian value (Equation 6.50)
        double J0 = A0(0,0) * A0(1,1) - A0(0,1) * A0(1,0);

        //Saint-Venant_Kirchhoff Stress
        //Equation 6.62 to 6.65)
        MatrixDouble S = fYoungModulus * E;

        //Terms from the second derivate of matrix A¹ (d²A¹/dtheta_beta dtheta_z)
        for(int alpha = 0; alpha < 3; alpha++){
            for(int beta = 0; beta < nphi; beta++){
                MatrixDouble DA1DY_beta(2,2);
                //Equation 6.69, 6.71 and 6.72
                ComputeDA1_DY(data, alpha, theta1, coordEta, beta, A0inv, A1, DA1DY_beta);
                
                //6.70 or 6.67
                MatrixDouble DEDy;
                ComputeDE_DY(A0inv, A1, DA1DY_beta, DEDy);

                int i = 3 * (beta) + alpha;
                double ScDEDY = DoubleContraction(S,DEDy);

                Rhs[i] -= (ScDEDY + forcingF[alpha]) * weightEta * J0 * data.fWeight * fDepth;
            }
        } 
    }
};

void PositionalFrame2D::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    PanicButton();
}

int PositionalFrame2D::VariableIndex(const std::string &name) const{
    
    if(!strcmp("Displacement",name.c_str()))           return 1;
    if(!strcmp("Rotation",name.c_str()))               return 2;
    if(!strcmp("ExactDisplacement",name.c_str()))      return 3;
    if(!strcmp("ExactRotation",name.c_str()))          return 4;
    if(!strcmp("ShearForce",name.c_str()))             return 5;
    if(!strcmp("BendingMoment",name.c_str()))          return 6;
    if(!strcmp("Load",name.c_str()))                   return 7;

    std::cout << "Post Process variable not implemented \n";
    PanicButton();
    return -1;
};

int PositionalFrame2D::NSolutionVariables(int var) const{
    switch (var)
    {
    case 1:
        return 3;
    case 3:
    case 5:
    case 2:
    case 4:
    case 6:
    case 7:
        return 1;

    default:
        PanicButton();
        return -1;
    }
};

void PositionalFrame2D::Solution(IntPointData &data, int var, VecDouble &Sol) {

    auto force = fForceFunction;
    VecDouble forcingF(3);
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    VecDouble phi,dphi,d2phi,d3phi;
    
    //Displacement
    if (var == 1){
        Sol[0] = data.fSol[0];
        Sol[1] = data.fSol[1];
        return;
    };

    //Rotation
    if (var == 2){
        Sol[0] = data.fSol[2];
        return;
    };

}; 

void PositionalFrame2D::ComputeA0A1(IntPointData &data, double coordEta, double &theta0, double &theta1, MatrixDouble &A0, MatrixDouble &A1){
    double dTheta0Dxi = 0.;

    //TODO: the following code works to straigt bars. Please implement the computation of theta_0 in the beggining of the analysis to ensure curved bars will be properly computed.
    for (int i = 0; i < data.fPhi.size(); i++){
        dTheta0Dxi += data.fDPhi(0,i) * theta0;
    }

    //Initial configuration mapping gradient (Equation 6.35 to 6.38)
    A0(0,0) = data.fAxes0(0,0) - 0.5 * fHeight * coordEta * sin(theta0) * dTheta0Dxi;
    A0(0,1) = 0.5 * fHeight * cos(theta0);
    A0(1,0) = data.fAxes0(1,0) + 0.5 * fHeight * coordEta * cos(theta0) * dTheta0Dxi;
    A0(1,1) = 0.5 * fHeight * sin(theta0);

    //Current configuration mapping gradient (Equation 6.39 to 6.42)
    A1(0,0) = data.fAxes1(0,0) - 0.5 * fHeight * coordEta * sin(theta1) * data.fDSolDAdim(2,0);
    A1(0,1) = 0.5 * fHeight * cos(theta1);
    A1(1,0) = data.fAxes1(1,0) + 0.5 * fHeight * coordEta * cos(theta1) * data.fDSolDAdim(2,0);
    A1(1,1) = 0.5 * fHeight * sin(theta1);

    data.fJacA0 = A0.determinant();

}

void PositionalFrame2D::ComputeGreenStrain(IntPointData &data, MatrixDouble &A0inv, MatrixDouble &A1, MatrixDouble &E){
    //Cauchy_Green Stretching (Equation 6.46)
    MatrixDouble A = A1 * A0inv;
    double detA = A.determinant();    
    //Green Deformation
    MatrixDouble Identity = Matrix2d::Identity();
    E = 0.5 * (A.transpose()*A - Identity);
}

void PositionalFrame2D::ComputeDA1_DY(IntPointData &data, int &direction, double &theta1, double &coordEta, int &beta, MatrixDouble &A0inv, MatrixDouble &A1, MatrixDouble &DA1_DY){
    DA1_DY.setZero();
    switch (direction){
    case 0:
        DA1_DY(0,0) = data.fDPhi(0, beta);
        break;
    case 1:
        DA1_DY(1,0) = data.fDPhi(0, beta);
        break;
    case 2:
        DA1_DY(0,0) = -0.5 * fHeight * coordEta * 
                            cos(theta1) * data.fPhi[beta] * data.fDSolDAdim(2,0) 
                            -0.5 * fHeight * coordEta *
                            sin(theta1) * data.fDPhi(0, beta);
        DA1_DY(1,0) = -0.5 * fHeight * coordEta * 
                            sin(theta1) * data.fPhi[beta] * data.fDSolDAdim(2,0) 
                            +0.5 * fHeight * coordEta *
                            cos(theta1) * data.fDPhi(0, beta);
        DA1_DY(0,1) = -0.5 * fHeight * (sin(theta1) * data.fPhi[beta]); 
        DA1_DY(1,1) =  0.5 * fHeight * (cos(theta1) * data.fPhi[beta]);
        break;
    
    default:
        PanicButton();
        break;
    }
}

void PositionalFrame2D::ComputeDE_DY(MatrixDouble &A0inv, MatrixDouble &A1, MatrixDouble &DA1_DY, MatrixDouble &DE_DY){
    DE_DY = 0.5 * (A0inv.transpose() * DA1_DY.transpose() * A1 * A0inv 
                 + A0inv.transpose() * A1.transpose() * DA1_DY * A0inv);
}


double PositionalFrame2D::DoubleContraction(const MatrixDouble &A, const MatrixDouble &B){
    return A(0,0)*B(0,0) + A(1,0)*B(1,0) + A(0,1)*B(0,1) + A(1,1)*B(1,1);
}