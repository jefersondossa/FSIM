#include "ElNavierStokes.h"

void ElNavierStokes::ComputeStiffness(int &index, MatrixDouble &Stiffness){
    double &visc_ = Mesh()->getProblemParameters().GetViscosity();
    double &dens_ = Mesh()->getProblemParameters().GetDensity();
    int DIM = Mesh()->Dimension();

    this->GetStabilizationParameter(index, this->tSUPG_, this->tPSPG_, this->tLSIC_, fIntegData.fDPhiX0);
    double WJ = fIntegData.fWeight * fIntegData.fJacA0 * getIntegPointWeightFunction(index);
    this->tSUPG_ = 0.;
    this->tLSIC_ = 0.;

    VecDouble u_(DIM+1);
    interpolateSolution(index, u_);

    //Solution Derivatives
    MatrixDouble du_dx(DIM+1,DIM);
    interpolateSolDerivatives(du_dx);
    auto dphi_dx = fIntegData.fDPhiX0;

    for (int i = Mesh()->NElNodes(); i-- ; ){        
        double shapeFi = Mesh()->getNumericalIntegration()-> phi_(i,index);
        double wSUPGi = 0.;
        for (int m=DIM; m--; ) wSUPGi += u_[m] * dphi_dx(i,m);
 
        for (int j = Mesh()->NElNodes(); j-- ; ){           

            double shapeFj = Mesh()->getNumericalIntegration()-> phi_(j,index);
            double shapeFij = shapeFi * shapeFj;
            double wSUPGj = 0.;
            for (int m=DIM; m--; ) wSUPGj += u_[m] * dphi_dx(j,m);
            
            //Convection matrix
            double C = (wSUPGj * shapeFi + wSUPGi * wSUPGj * tSUPG_) * dens_ * 0.;

            double aux1 = tSUPG_ * wSUPGi * shapeFj;
            double aux2 = tSUPG_ * shapeFj;

            for (int k = DIM; k--;  ){
                Stiffness((DIM+1)*i+k,(DIM+1)*j+k) += C * WJ;
                double conv = 0.;
                for (int m = DIM; m--; ) conv += u_[m]*du_dx(k,m);

                for (int l = DIM; l--; ){

                    //Diffusion matrix
                    double K = fIntegData.fDPhiX0(i,l) * fIntegData.fDPhiX0(j,k) * visc_;
                    if (k==l) for (int m = DIM; m--; ) K += fIntegData.fDPhiX0(i,m) * fIntegData.fDPhiX0(j,m)* visc_;

                    //Convection derivatives
                    double Cuu = (shapeFij * du_dx(k,l) + 
                                  aux1 * du_dx(k,l) +
                                  aux2 * conv * fIntegData.fDPhiX0(i,l)) * dens_;

                    //LSIC
                    double KLS = fIntegData.fDPhiX0(i,k) * fIntegData.fDPhiX0(j,l) * tLSIC_ * dens_;

                    Stiffness((DIM+1)*i+k,(DIM+1)*j+l) += (K + KLS + Cuu) * WJ;
                }
                //Gradient operator
                double Q_SUPG = - fIntegData.fDPhiX0(i,k) * shapeFj + wSUPGi * fIntegData.fDPhiX0(j,k) * tSUPG_;
                //Divergent operator
                double Q = fIntegData.fDPhiX0(i,k) * shapeFj;

                Stiffness((DIM+1)*i+k,(DIM+1)*j+DIM) += Q_SUPG * WJ;
                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += Q * WJ;

                //PSPG stabilization
                double G = fIntegData.fDPhiX0(i,k) * wSUPGj * tPSPG_;
                double Guu = 0.;
                for (int m = DIM; m--; ) Guu += fIntegData.fDPhiX0(i,m) * du_dx(m,k) * shapeFj * tPSPG_;

                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += (G + Guu) * WJ;
            }
            //PSPG stabilization
            double Q = 0.;
            for (int m = DIM; m--; ) Q += fIntegData.fDPhiX0(i,m) * fIntegData.fDPhiX0(j,m) * tPSPG_ / dens_;
            Stiffness((DIM+1)*j+DIM,(DIM+1)*i+DIM) += Q * WJ;
        };
    };



    for (int i = Mesh()->NElNodes(); i-- ; ){        
        double shapeFi = Mesh()->getNumericalIntegration()-> phi_(i,index);
        for (int j = Mesh()->NElNodes(); j-- ; ){
            
            double shapeFj = Mesh()->getNumericalIntegration()-> phi_(j,index);

            for (int k = DIM; k--;  ){
                for (int l = DIM; l--; ){

                    //Diffusion matrix
                    double K = fIntegData.fDPhiX0(i,l) * fIntegData.fDPhiX0(j,k) * visc_;
                    if (k==l) for (int m = DIM; m--; ) K += fIntegData.fDPhiX0(i,m) * fIntegData.fDPhiX0(j,m)* visc_;

                    Stiffness((DIM+1)*i+k,(DIM+1)*j+l) += K * WJ;
                }
                //Gradient operator
                //Divergent operator
                double Q = fIntegData.fDPhiX0(i,k) * shapeFj;

                Stiffness((DIM+1)*i+k,(DIM+1)*j+DIM) += -Q * WJ;
                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += Q * WJ;
            }
            //PSPG stabilization
            double Q = 0.;
            for (int m = DIM; m--; ) Q += fIntegData.fDPhiX0(i,m) * fIntegData.fDPhiX0(j,m) * tPSPG_ / dens_;
            Stiffness((DIM+1)*j+DIM,(DIM+1)*i+DIM) += Q * WJ;
        };
    };

}

void ElNavierStokes::ComputeResidual(int &index, VecDouble &Rhs){
    VecDouble fieldForce = Mesh()->getProblemParameters().GetFieldForce();
    auto force = Mesh()->getProblemParameters().getForcingFunction();
    int DIM = Mesh()->Dimension();
    double &visc_ = Mesh()->getProblemParameters().GetViscosity();
    double &dens_ = Mesh()->getProblemParameters().GetDensity();

    VecDouble forcingF(DIM);
    forcingF.setZero();
    VecDouble x_ = getIntegPointCoordinatesValue(index);
    if (force) force(x_,forcingF);

    //Solution Derivatives
    MatrixDouble du_dx(DIM+1,DIM);
    interpolateSolDerivatives(du_dx);

    VecDouble u_(DIM+1);
    interpolateSolution(index, u_);

    double WJ = fIntegData.fWeight * fIntegData.fJacA0 * getIntegPointWeightFunction(index);

    double divrU = 0.;
    for (int l=DIM; l--; ) divrU += du_dx(l,l);

    for (int i = Mesh()->NElNodes(); i--; ){
        // std::cout << "Sol = " << Mesh()->NodeVec()[getConnectivity()[i]]->GetSolution(0) << std::endl;
        double shapeFi = Mesh()->getNumericalIntegration()-> phi_(i,index);

        for (int k = DIM; k--; ){
            //Viscosity
            double K = 0.;
            for (int l=DIM; l--; ) K += fIntegData.fDPhiX0(i,l) * du_dx(k,l) * visc_;
            for (int l=DIM; l--; ) K += fIntegData.fDPhiX0(i,l) * du_dx(l,k) * visc_;

            //LSIC
            double KLS = fIntegData.fDPhiX0(i,k) * divrU * tLSIC_ * dens_;

            //Convection + SUPG
            double C = 0.;
            for (int l=DIM; l--; ) C += du_dx(k,l) * u_[l] * shapeFi * dens_;
            double conv = 0.;
            for (int l=DIM; l--; ) conv += u_[l] * fIntegData.fDPhiX0(i,l);
            for (int l=DIM; l--; ) C += conv * u_[l] * du_dx(k,l) * tSUPG_ * dens_;

            //Pressure + SUPG
            double P = - (fIntegData.fDPhiX0(i,k) * u_[DIM]);
            for (int l=DIM; l--; ) P += fIntegData.fDPhiX0(i,l) * u_[l] * du_dx(DIM,k) * tSUPG_;

            //External force
            double F = (fieldForce[k]*dens_ + forcingF[k]) * shapeFi;
            
            Rhs[(DIM+1)*i+k] += (-K - P - C - KLS + F) * WJ;
        }

        double Q = divrU * shapeFi;
        for (int l=DIM; l--; ) Q += fIntegData.fDPhiX0(i,l) * du_dx(DIM,l) * tPSPG_ / dens_
                                  + fIntegData.fDPhiX0(i,l) * (fieldForce[l] + forcingF[l]/dens_) * tPSPG_;
        for (int k=DIM; k--; )
            for (int l=DIM; l--; )
                Q += fIntegData.fDPhiX0(i,k) * u_[l] * du_dx(k,l) * tPSPG_;

        Rhs[(DIM+1)*i+DIM] += -Q * WJ;
                            
    };
    
};

void ElNavierStokes::ComputeError(VecDouble &errors){
    std::cout << "Not implemented yet.\n";
    // PanicButton();
}
