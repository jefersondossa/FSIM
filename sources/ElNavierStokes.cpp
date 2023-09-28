#include "ElNavierStokes.h"

void ElNavierStokes::ComputeStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, MatrixDouble &Stiffness){
    double &visc_ = Mesh()->getProblemParameters().GetViscosity();
    double &dens_ = Mesh()->getProblemParameters().GetDensity();
    int DIM = Mesh()->Dimension();

    this->GetStabilizationParameter(index, this->tSUPG_, this->tPSPG_, this->tLSIC_, dphi_dx);
    double WJ = weight_ * djac_ * getIntegPointWeightFunction(index);

    VecDouble u_(DIM+1);
    interpolateSolution(index, u_);

    //Solution Derivatives
    MatrixDouble du_dx(DIM+1,DIM);
    interpolateSolDerivatives(dphi_dx, du_dx);

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
            double C = (wSUPGj * shapeFi + wSUPGi * wSUPGj * tSUPG_) * dens_;

            double aux1 = tSUPG_ * wSUPGi * shapeFj;
            double aux2 = tSUPG_ * shapeFj;

            for (int k = DIM; k--;  ){
                Stiffness((DIM+1)*i+k,(DIM+1)*j+k) += C * WJ;
                double conv = 0.;
                for (int m = DIM; m--; ) conv += u_[m]*du_dx(k,m);

                for (int l = DIM; l--; ){

                    //Diffusion matrix
                    double K = dphi_dx(i,l) * dphi_dx(j,k) * visc_;
                    if (k==l) for (int m = DIM; m--; ) K += dphi_dx(i,m) * dphi_dx(j,m)* visc_;

                    //Convection derivatives
                    double Cuu = (shapeFij * du_dx(k,l) + 
                                  aux1 * du_dx(k,l) +
                                  aux2 * conv * dphi_dx(i,l)) * dens_*0.;

                    //LSIC
                    double KLS = dphi_dx(i,k) * dphi_dx(j,l) * tLSIC_ * dens_;

                    Stiffness((DIM+1)*i+k,(DIM+1)*j+l) += (K + KLS + Cuu) * WJ;
                }
                //Gradient operator
                double Q_SUPG = - dphi_dx(i,k) * shapeFj + wSUPGi * dphi_dx(j,k) * tSUPG_;
                //Divergent operator
                double Q = dphi_dx(i,k) * shapeFj;

                Stiffness((DIM+1)*i+k,(DIM+1)*j+DIM) += Q_SUPG * WJ;
                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += Q * WJ;

                //PSPG stabilization
                double G = dphi_dx(i,k) * wSUPGj * tPSPG_;
                double Guu = 0.;
                for (int m = DIM; m--; ) Guu += dphi_dx(i,m) * du_dx(m,k) * shapeFj * tPSPG_*0.;

                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += (G + Guu) * WJ;
            }
            //PSPG stabilization
            double Q = 0.;
            for (int m = DIM; m--; ) Q += dphi_dx(i,m) * dphi_dx(j,m) * tPSPG_ / dens_;
            Stiffness((DIM+1)*j+DIM,(DIM+1)*i+DIM) += Q * WJ;
        };
    };

}

void ElNavierStokes::ComputeResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, VecDouble &Rhs){
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
    interpolateSolDerivatives(dphi_dx, du_dx);

    VecDouble u_(DIM+1);
    interpolateSolution(index, u_);

    double WJ = weight_ * djac_ * getIntegPointWeightFunction(index);

    double divrU = 0.;
    for (int l=DIM; l--; ) divrU += du_dx(l,l);

    for (int i = Mesh()->NElNodes(); i--; ){
        // std::cout << "Sol = " << Mesh()->NodeVec()[getConnectivity()[i]]->GetSolution(0) << std::endl;
        double shapeFi = Mesh()->getNumericalIntegration()-> phi_(i,index);

        for (int k = DIM; k--; ){
            //Viscosity
            double K = 0.;
            for (int l=DIM; l--; ) K += dphi_dx(i,l) * du_dx(k,l) * visc_;
            for (int l=DIM; l--; ) K += dphi_dx(i,l) * du_dx(l,k) * visc_;

            //LSIC
            double KLS = dphi_dx(i,k) * divrU * tLSIC_ * dens_;

            //Convection + SUPG
            double C = 0.;
            for (int l=DIM; l--; ) C += du_dx(k,l) * u_[l] * shapeFi * dens_;
            double conv = 0.;
            for (int l=DIM; l--; ) conv += u_[l] * dphi_dx(i,l);
            for (int l=DIM; l--; ) C += conv * u_[l] * du_dx(k,l) * tSUPG_ * dens_*0.;

            //Pressure + SUPG
            double P = - (dphi_dx(i,k) * u_[DIM]);
            for (int l=DIM; l--; ) P += dphi_dx(i,l) * u_[l] * du_dx(DIM,k) * tSUPG_;

            //External force
            double F = (fieldForce[k]*dens_ + forcingF[k]) * shapeFi;
            
            Rhs[(DIM+1)*i+k] += (-K - P - C - KLS + F) * WJ;
        }

        double Q = divrU * shapeFi;
        for (int l=DIM; l--; ) Q += dphi_dx(i,l) * du_dx(DIM,l) * tPSPG_ / dens_
                                  + dphi_dx(i,l) * (fieldForce[l] + forcingF[l]/dens_) * tPSPG_;
        for (int k=DIM; k--; )
            for (int l=DIM; l--; )
                Q += dphi_dx(i,k) * u_[l] * du_dx(k,l) * tPSPG_*0.;

        Rhs[(DIM+1)*i+DIM] += -Q * WJ;
                            
    };
    
};

void ElNavierStokes::ComputeError(VecDouble &errors){
    std::cout << "Not implemented yet.\n";
    // PanicButton();
}
