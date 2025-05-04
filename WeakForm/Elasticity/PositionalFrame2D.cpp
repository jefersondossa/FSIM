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

    IntRule1d ruleEta(6,3);

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

    double dTheta0Dxi = 0.;

    //TODO: the following code works to straigt bars. Please implement the computation of theta_0 in the beggining of the analysis to ensure curved bars will be properly computed.
    for (int i = 0; i < nphi; i++){
        dTheta0Dxi += data.fDPhi(0,i) * theta0;
    }
    
    for (int eta = 0; eta < ruleEta.NPoints(); eta++){
        double coordEta = ruleEta.PointList(eta,0);
        double weightEta = ruleEta.WeightList(eta);
        

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

        MatrixDouble A0inv = A0.inverse();

        //Cauchy_Green Stretching (Equation 6.46)
        MatrixDouble C = A0inv.transpose() * A1.transpose() * A1 * A0inv;
        
        //Green Deformation
        MatrixDouble E(2,2);
        MatrixDouble Identity = Matrix2d::Identity();

        E = 0.5 * (C - Identity);

        //Jacobian value (Equation 6.50)
        double J0 = A0(0,0) * A0(1,1) - A0(0,1) * A0(1,0);

        //Saint-Venant_Kirchhoff Stress
        MatrixDouble S(2,2);

        //Equation 6.62 to 6.65)
        S = fYoungModulus * E;
        
        // Variable initialization
        std::vector<std::vector<MatrixDouble>> DA1DYbeta(3), DEDy(3);

        DA1DYbeta[0].resize(nphi);
        DA1DYbeta[1].resize(nphi);
        DA1DYbeta[2].resize(nphi);

        DEDy[0].resize(nphi);
        DEDy[1].resize(nphi);
        DEDy[2].resize(nphi);

        for (int beta = 0; beta < nphi; beta++){

            //Equation 6.69, 6.71 and 6.72
            DA1DYbeta[0][beta].resize(2,2);
            DA1DYbeta[1][beta].resize(2,2);
            DA1DYbeta[2][beta].resize(2,2);
            
            DA1DYbeta[0][beta].setZero();
            DA1DYbeta[1][beta].setZero();
            DA1DYbeta[2][beta].setZero();

            DA1DYbeta[0][beta](0,0) = data.fDPhi(0, beta);
            DA1DYbeta[1][beta](1,0) = data.fDPhi(0, beta);
            
            DA1DYbeta[2][beta](0,0) = 0.5 * fHeight * coordEta *
                                     (-cos(theta1) * data.fPhi[beta] * data.fDSolDAdim(2,0) 
                                      -sin(theta1) * data.fDPhi(0, beta));
            DA1DYbeta[2][beta](1,0) = 0.5 * fHeight * coordEta * 
                                     (-sin(theta1) * data.fPhi[beta] * data.fDSolDAdim(2,0) 
                                      +cos(theta1) * data.fDPhi(0, beta));
            DA1DYbeta[2][beta](0,1) = -0.5 * fHeight * (sin(theta1) * data.fPhi[beta]); 
            DA1DYbeta[2][beta](1,1) =  0.5 * fHeight * (cos(theta1) * data.fPhi[beta]);
            //6.70 or 6.67
            DEDy[0][beta] = 0.5 * (A0inv.transpose() * DA1DYbeta[0][beta].transpose() * A1 * A0inv 
                                 + A0inv.transpose() *  A1.transpose() *DA1DYbeta[0][beta] * A0inv);
            DEDy[1][beta] = 0.5 * (A0inv.transpose() * DA1DYbeta[1][beta].transpose() * A1 * A0inv 
                                 + A0inv.transpose() *  A1.transpose() *DA1DYbeta[1][beta] * A0inv);
            DEDy[2][beta] = 0.5 * (A0inv.transpose() * DA1DYbeta[2][beta].transpose() * A1 * A0inv 
                                 + A0inv.transpose() *  A1.transpose() *DA1DYbeta[2][beta] * A0inv);  
            
            // PrintMathematica(DEDy[2][beta]);
        }

        //Terms from the second derivate of matrix A¹ (d²A¹/dtheta_beta dtheta_z)
        MatrixDouble D2DA1(2,2);
        for(int alpha = 0; alpha < 3; alpha++){
            for(int beta = 0; beta < nphi; beta++){
                for (int gama = 0; gama < 3; gama++){
                    for (int zeta = 0; zeta < nphi; zeta++){
                        
                        //Equation 6.83 to 6.86
                        D2DA1(0,0) = 0.5 * fHeight * coordEta * (sin(theta1) * data.fPhi[zeta] * data.fPhi[beta] * data.fDSolDAdim(2,0) 
                                                               - cos(theta1) * (data.fPhi[zeta] * data.fDPhi(0,beta) + data.fPhi[beta] * data.fDPhi(0,zeta)));
                        D2DA1(1,0) =-0.5 * fHeight * coordEta * (cos(theta1) * data.fPhi[zeta] * data.fPhi[beta] * data.fDSolDAdim(2,0) 
                                                               + sin(theta1) * (data.fPhi[zeta] * data.fDPhi(0,beta) + data.fPhi[beta] * data.fDPhi(0,zeta)));
                        D2DA1(0,1) = -0.5 * fHeight * (cos(theta1) * data.fPhi[zeta] * data.fPhi[beta]); 
                        D2DA1(1,1) = -0.5 * fHeight * (sin(theta1) * data.fPhi[zeta] * data.fPhi[beta]);

                        //Equation 6.81
                        MatrixDouble temp1 = A0inv.transpose() * DA1DYbeta[alpha][beta].transpose() * DA1DYbeta[gama][zeta]* A0inv;
                        MatrixDouble temp2 = A0inv.transpose() * D2DA1.transpose() * A1 * A0inv;
                        MatrixDouble D2EDY2 = 0.5 * (temp1 + temp1.transpose() + temp2 + temp2.transpose());
                        // MatrixDouble D2EDY2 = 0.5 * (A0inv.transpose() * DA1DYbeta[alpha][beta].transpose() * DA1DYbeta[gama][zeta]* A0inv + 
                        //                              A0inv.transpose() * DA1DYbeta[gama][zeta].transpose() * DA1DYbeta[alpha][beta]* A0inv) +
                        //                             (A0inv.transpose() * D2DA1.transpose() * A1 * A0inv + 
                        //                              A0inv.transpose() * A1.transpose() * A1 * A0inv);
                                                    //  !Calculating the 2nd derivative of the Green strain tensor:

                                                    //  temp1 = matmul( matmul(transpose(A0_inv) , transpose(self%calc_1der_A1(nd1,dof1,ig_rl,ig_h))) , matmul(self%calc_1der_A1(nd2,dof2,ig_rl,ig_h) , A0_inv) )
                                                    //  temp2 = matmul( matmul(transpose(A0_inv) , transpose(d2A1_dY2)) , matmul(self%kin_intpt(ig_rl,ig_h)%A1 , A0_inv) )
                                                    //  calc_2der_Green_strain = 0.5d0 * (temp1 + transpose(temp1) + temp2 + transpose(temp2))
                        //Equation 6.79
                        MatrixDouble DSDy = fYoungModulus * DEDy[gama][zeta];

                        //Equation 6.76
                        double DSDYcDEDy = DSDy(0,0)*DEDy[alpha][beta](0,0) 
                                         + DSDy(1,0)*DEDy[alpha][beta](1,0) 
                                         + DSDy(0,1)*DEDy[alpha][beta](0,1) 
                                         + DSDy(1,1)*DEDy[alpha][beta](1,1);
                        double ScD2EDY2 = S(0,0) * D2EDY2(0,0) 
                                        + S(1,0) * D2EDY2(1,0) 
                                        + S(0,1) * D2EDY2(0,1) 
                                        + S(1,1) * D2EDY2(1,1);                        

                        int j = 3 * (beta) + alpha;
                        int i = 3 * (zeta) + gama;
                        // std::cout << "i = " << i << ", j = " << j << std::endl;
                        // std::cout << "DSDYcDEDy = " << DSDYcDEDy << ", ScD2EDY2 = " << ScD2EDY2 << std::endl;

                        Stiffness(i,j) += (DSDYcDEDy+ScD2EDY2) * fDepth * weightEta * J0 * data.fWeight;
                        // std::cout << "pesos = " << weightEta*data.fWeight << std::endl;
                        
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

    IntRule1d ruleEta(6,3);

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

    double dTheta0Dxi = 0.;

    //TODO: the following code works to straigt bars. Please implement the computation of theta_0 in the beggining of the analysis to ensure curved bars will be properly computed.
    for (int i = 0; i < nphi; i++){
        dTheta0Dxi += data.fDPhi(0,i) * theta0;
    }
    
    for (int eta = 0; eta < ruleEta.NPoints(); eta++){
        double coordEta = ruleEta.PointList(eta,0);
        double weightEta = ruleEta.WeightList(eta);
        

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

        MatrixDouble A0inv = A0.inverse();

        //Cauchy_Green Stretching (Equation 6.46)
        MatrixDouble C = A0inv.transpose() * A1.transpose() * A1 * A0inv;
        
        //Green Deformation
        MatrixDouble E(2,2);
        MatrixDouble Identity = Matrix2d::Identity();

        E = 0.5 * (C - Identity);

        //Jacobian value (Equation 6.50)
        double J0 = A0(0,0) * A0(1,1) - A0(0,1) * A0(1,0);

        //Saint-Venant_Kirchhoff Stress
        //Equation 6.62 to 6.65)
        MatrixDouble S = fYoungModulus * E;

        //Terms from the second derivate of matrix A¹ (d²A¹/dtheta_beta dtheta_z)
        for(int alpha = 0; alpha < 3; alpha++){
            for(int beta = 0; beta < nphi; beta++){
                MatrixDouble DA1DYbeta(2,2);
                DA1DYbeta.setZero();
                //Equation 6.69, 6.71 and 6.72
                switch (alpha){
                case 0:
                    DA1DYbeta(0,0) = data.fDPhi(0, beta);
                    break;
                case 1:
                    DA1DYbeta(1,0) = data.fDPhi(0, beta);
                    break;
                case 2:
                    DA1DYbeta(0,0) = 0.5 * fHeight * coordEta *
                                    (-cos(theta1) * data.fPhi[beta] * data.fDSolDAdim(2,0) 
                                    -sin(theta1) * data.fDPhi(0, beta));
                    DA1DYbeta(1,0) = 0.5 * fHeight * coordEta * 
                                    (-sin(theta1) * data.fPhi[beta] * data.fDSolDAdim(2,0) 
                                    +cos(theta1) * data.fDPhi(0, beta));
                    DA1DYbeta(0,1) = -0.5 * fHeight * (sin(theta1) * data.fPhi[beta]); 
                    DA1DYbeta(1,1) =  0.5 * fHeight * (cos(theta1) * data.fPhi[beta]);
                    break;
                
                default:
                    PanicButton();
                    break;
                }
                //6.70 or 6.67
                MatrixDouble DEDy = 0.5 * 
                            (A0inv.transpose() * DA1DYbeta.transpose() * A1 * A0inv 
                           + A0inv.transpose() *  A1.transpose() *DA1DYbeta * A0inv);
                
                int i = 3 * (beta) + alpha;
                // std::cout << "i = " << i << ", j = " << j << std::endl;
                // std::cout << "DSDYcDEDy = " << DSDYcDEDy << ", ScD2EDY2 = " << ScD2EDY2 << std::endl;
                double ScDEDY = S(0,0) * DEDy(0,0) 
                              + S(1,0) * DEDy(1,0) 
                              + S(0,1) * DEDy(0,1) 
                              + S(1,1) * DEDy(1,1);
                // if (i==5){
                //     std::cout << "i = " << i << ", alpha = " << alpha << " , beta = " << beta << std::endl;
                //     PrintMathematica(S);
                //     PrintMathematica(DEDy[alpha][beta]);
                //     std::cout << "DSDYcDEDy = " << ScDEDY << std::endl;                   
                //     std::cout << weightEta*data.fWeight*J0 << std::endl;
                // }            

                Rhs[i] += (-ScDEDY - forcingF[alpha]) * weightEta * J0 * data.fWeight * fDepth;
                
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

    HermiteFunction(data.fAdimCoord[0],phi,dphi,d2phi,d3phi);
    
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

void PositionalFrame2D::HermiteFunction(double &ksi, VecDouble &phi, VecDouble &dphi, VecDouble &d2phi, VecDouble &d3phi){
    phi.resize(4);
    phi[0] = 0.5 - 0.75*ksi + 0.25*ksi*ksi*ksi;
    phi[1] = 0.25 - 0.25*ksi - 0.25*ksi*ksi + 0.25*ksi*ksi*ksi;
    phi[2] = 0.5 + 0.75*ksi - 0.25*ksi*ksi*ksi;
    phi[3] = -0.25 - 0.25*ksi + 0.25*ksi*ksi + 0.25*ksi*ksi*ksi;
    dphi.resize(4);
    dphi[0] = -0.75 + 0.75*ksi*ksi;
    dphi[1] = -0.25 - 0.5*ksi + 0.75*ksi*ksi;
    dphi[2] = 0.75 - 0.75*ksi*ksi;
    dphi[3] = -0.25 + 0.5*ksi + 0.75*ksi*ksi;
    d2phi.resize(4);
    d2phi[0] = 1.5*ksi;
    d2phi[1] = -0.5 + 1.5*ksi;
    d2phi[2] = -1.5*ksi;
    d2phi[3] = 0.5 + 1.5*ksi;
    d3phi.resize(4);
    d3phi[0] = 1.5;
    d3phi[1] = 1.5;
    d3phi[2] = -1.5;
    d3phi[3] = 1.5;

}