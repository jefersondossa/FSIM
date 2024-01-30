#include "NavierStokes.h"

NavierStokes::NavierStokes(int matid, int dim, double density, double viscosity) : Stokes(matid,dim,density,viscosity) {
    
};

void NavierStokes::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    
    if (!data.fNeedsDSol || !data.fNeedsSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
        data.fNeedsSol = true;
        data.fSol.resize(fNState);
    }
    int DIM = fDimension;
    auto dphi_dx = data.fDPhiX0;
    this->GetStabilizationParameter(index, data);
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    // this->tSUPG_ = 0.;
    // this->tLSIC_ = 0.;
    int nphi = data.fPhi.size();

    for (int i = nphi; i-- ; ){        
        double shapeFi = data.fPhi[i];
        double wSUPGi = 0.;
        for (int m=DIM; m--; ) wSUPGi += data.fSol[m] * dphi_dx(i,m);
 
        for (int j = nphi; j-- ; ){           

            double shapeFj = data.fPhi[j];
            double shapeFij = shapeFi * shapeFj;
            double wSUPGj = 0.;
            for (int m=DIM; m--; ) wSUPGj += data.fSol[m] * dphi_dx(j,m);
            
            //Convection matrix
            double C = (wSUPGj * shapeFi + wSUPGi * wSUPGj * this->tSUPG_) * fDensity * 0.;

            double aux1 = this->tSUPG_ * wSUPGi * shapeFj;
            double aux2 = this->tSUPG_ * shapeFj;

            for (int k = DIM; k--;  ){
                Stiffness((DIM+1)*i+k,(DIM+1)*j+k) += C * WJ;
                double conv = 0.;
                for (int m = DIM; m--; ) conv += data.fSol[m]*data.fDSolDx(k,m);

                for (int l = DIM; l--; ){

                    //Diffusion matrix
                    double K = data.fDPhiX0(i,l) * data.fDPhiX0(j,k) * fViscosity;
                    if (k==l) for (int m = DIM; m--; ) K += data.fDPhiX0(i,m) * data.fDPhiX0(j,m)* fViscosity;

                    //Convection derivatives
                    double Cuu = (shapeFij * data.fDSolDx(k,l) + 
                                  aux1 * data.fDSolDx(k,l) +
                                  aux2 * conv * data.fDPhiX0(i,l)) * fDensity;

                    //LSIC
                    double KLS = data.fDPhiX0(i,k) * data.fDPhiX0(j,l) * this->tLSIC_ * fDensity;

                    Stiffness((DIM+1)*i+k,(DIM+1)*j+l) += (K + KLS + Cuu) * WJ;
                }
                //Gradient operator
                double Q_SUPG = - data.fDPhiX0(i,k) * shapeFj + wSUPGi * data.fDPhiX0(j,k) * this->tSUPG_;
                //Divergent operator
                double Q = data.fDPhiX0(i,k) * shapeFj;

                Stiffness((DIM+1)*i+k,(DIM+1)*j+DIM) += Q_SUPG * WJ;
                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += Q * WJ;

                //PSPG stabilization
                double G = data.fDPhiX0(i,k) * wSUPGj * this->tPSPG_;
                double Guu = 0.;
                for (int m = DIM; m--; ) Guu += data.fDPhiX0(i,m) * data.fDSolDx(m,k) * shapeFj * this->tPSPG_;

                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += (G + Guu) * WJ;
            }
            //PSPG stabilization
            double Q = 0.;
            for (int m = DIM; m--; ) Q += data.fDPhiX0(i,m) * data.fDPhiX0(j,m) * this->tPSPG_ / fDensity;
            Stiffness((DIM+1)*j+DIM,(DIM+1)*i+DIM) += Q * WJ;
        };
    };



    for (int i = nphi; i-- ; ){        
        double shapeFi = data.fPhi[i];
        for (int j = nphi; j-- ; ){
            
            double shapeFj = data.fPhi[j];

            for (int k = DIM; k--;  ){
                for (int l = DIM; l--; ){

                    //Diffusion matrix
                    double K = data.fDPhiX0(i,l) * data.fDPhiX0(j,k) * fViscosity;
                    if (k==l) for (int m = DIM; m--; ) K += data.fDPhiX0(i,m) * data.fDPhiX0(j,m)* fViscosity;

                    Stiffness((DIM+1)*i+k,(DIM+1)*j+l) += K * WJ;
                }
                //Gradient operator
                //Divergent operator
                double Q = data.fDPhiX0(i,k) * shapeFj;

                Stiffness((DIM+1)*i+k,(DIM+1)*j+DIM) += -Q * WJ;
                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += Q * WJ;
            }
            //PSPG stabilization
            double Q = 0.;
            for (int m = DIM; m--; ) Q += data.fDPhiX0(i,m) * data.fDPhiX0(j,m) * this->tPSPG_ / fDensity;
            Stiffness((DIM+1)*j+DIM,(DIM+1)*i+DIM) += Q * WJ;
        };
    };

}

void NavierStokes::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    auto force = fForceFunction;
    int DIM = fDimension;

    int nphi = data.fPhi.size();
    VecDouble forcingF(DIM);
    forcingF.setZero();
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];

    double divrU = 0.;
    for (int l=DIM; l--; ) divrU += data.fDSolDx(l,l);

    for (int i = nphi; i--; ){
        double shapeFi = data.fPhi[i];

        for (int k = DIM; k--; ){
            //Viscosity
            double K = 0.;
            for (int l=DIM; l--; ) K += data.fDPhiX0(i,l) * data.fDSolDx(k,l) * fViscosity;
            for (int l=DIM; l--; ) K += data.fDPhiX0(i,l) * data.fDSolDx(l,k) * fViscosity;

            //LSIC
            double KLS = data.fDPhiX0(i,k) * divrU * this->tLSIC_ * fDensity;

            //Convection + SUPG
            double C = 0.;
            for (int l=DIM; l--; ) C += data.fDSolDx(k,l) * data.fSol[l] * shapeFi * fDensity;
            double conv = 0.;
            for (int l=DIM; l--; ) conv += data.fSol[l] * data.fDPhiX0(i,l);
            for (int l=DIM; l--; ) C += conv * data.fSol[l] * data.fDSolDx(k,l) * this->tSUPG_ * fDensity;

            //Pressure + SUPG
            double P = - (data.fDPhiX0(i,k) * data.fSol[DIM]);
            for (int l=DIM; l--; ) P += data.fDPhiX0(i,l) * data.fSol[l] * data.fDSolDx(DIM,k) * this->tSUPG_;

            //External force
            double F = (forcingF[k]) * shapeFi;
            
            Rhs[(DIM+1)*i+k] += (-K - P - C - KLS + F) * WJ;
        }

        double Q = divrU * shapeFi;
        for (int l=DIM; l--; ) Q += data.fDPhiX0(i,l) * data.fDSolDx(DIM,l) * this->tPSPG_ / fDensity
                                  + data.fDPhiX0(i,l) * (forcingF[l]/fDensity) * this->tPSPG_;
        for (int k=DIM; k--; )
            for (int l=DIM; l--; )
                Q += data.fDPhiX0(i,k) * data.fSol[l] * data.fDSolDx(k,l) * this->tPSPG_;

        Rhs[(DIM+1)*i+DIM] += -Q * WJ;
                            
    };
    
};

void NavierStokes::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet.\n";
    // PanicButton();
}


int NavierStokes::VariableIndex(const std::string &name) const{
    PanicButton();
};

int NavierStokes::NSolutionVariables(int var) const{
    PanicButton();
};

void NavierStokes::Solution(IntPointData &data, int var, VecDouble &Sol){
    PanicButton();
};