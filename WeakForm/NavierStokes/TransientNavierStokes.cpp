#include "TransientNavierStokes.h"

TransientNavierStokes::TransientNavierStokes(int matid, int dim, double density, double viscosity, double dTime, double specRadius) : NavierStokes(matid,dim,density,viscosity) {
    fTimeStep = dTime;
    fIntegScheme = EGeneralizedAlpha;
    fSpectralRadius = specRadius;
    fAlphaF = 1. / (1. + specRadius);
    fAlphaM = 0.5 * (3. - specRadius)/(1. + specRadius);
    fGamma = 0.5 + fAlphaM - fAlphaF;
};

void TransientNavierStokes::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    //Velocity
    VecDouble una_ = fAlphaF * data.fSol + (1. - fAlphaF) * data.fSolPrev;

    //Velocity Derivatives
    MatrixDouble duna_dx = fAlphaF * data.fDSolDx + (1. - fAlphaF) * data.fDSolDxPrev;
    
    //Mesh Velocity
    VecDouble umesh_(fDimension), umeshPrev_(fDimension), umeshna_(fDimension);
    umeshna_ = fAlphaF * umesh_ + (1. - fAlphaF) * umeshPrev_;

    // Trust me, it improves performance!
    double AGDT = fAlphaF * fGamma * fTimeStep; 
    double VAGDT = fViscosity * AGDT; 
    double DAGDT = fDensity * AGDT; 
    double WJ = data.fWeight * data.fJacA0;
    double DAM = fDensity * fAlphaM;

    int nphi = data.fPhi.size();
    for (int i = nphi; i-- ; ){
        
        double shapeFi = data.fPhi[i];
        double wSUPGi = 0.;
        for (int m=fDimension; m--; ) wSUPGi += (una_[m] - umeshna_[m]) * data.fDPhiX0(i,m);

        for (int j = nphi; j-- ; ){
            
            double shapeFj = data.fPhi[j];
            double shapeFij = shapeFi * shapeFj;
            double wSUPGj = 0.;
            for (int m=fDimension; m--; ) wSUPGj += (una_[m] - umeshna_[m]) * data.fDPhiX0(j,m);
            
            //Mass matrix (use for both directions)
            double M = (shapeFij + wSUPGi * shapeFj * tSUPG_) * DAM;

            //Convection matrix
            double C = (wSUPGj * shapeFi + wSUPGi * wSUPGj * tSUPG_) * DAGDT;         

            double aux1 = tSUPG_ * wSUPGi * shapeFj;
            double aux2 = tSUPG_ * shapeFj;

            for (int k = fDimension; k--;  ){

                Stiffness(fDimension*i+k,fDimension*j+k) += (M + C) * WJ;
                    
                double conv = 0.;
                for (int m = fDimension; m--; ) conv += una_[m]*duna_dx(k,m);

                for (int l = fDimension; l--; ){

                    //Diffusion matrix
                    double K = data.fDPhiX0(i,l) * data.fDPhiX0(j,k) * VAGDT;
                    if (k==l) for (int m = fDimension; m--; ) K += data.fDPhiX0(i,m) * data.fDPhiX0(j,m) * VAGDT;

                    //Convection derivatives
                    double Cuu = (shapeFij * duna_dx(k,l) + 
                                  aux1 * duna_dx(k,l) +
                                  aux2 * conv * data.fDPhiX0(i,l)) * DAGDT;

                    //LSIC
                    double KLS = data.fDPhiX0(i,k) * data.fDPhiX0(j,l) * tLSIC_ * DAGDT;

                    Stiffness(fDimension*i+k,fDimension*j+l) += (K + KLS + Cuu) * WJ;
                }

                //SINAL DA PARCELA QUE MULTIPLICA O TSUPG ESTA COM SINAL TROCADO NA FORMULAÇAO DO TEZDUYAR
                //Gradient operator
                double Q_SUPG = - data.fDPhiX0(i,k) * shapeFj + wSUPGi * data.fDPhiX0(j,k) * tSUPG_;
                //Divergent operator
                double Q = data.fDPhiX0(i,k) * shapeFj * AGDT;

                Stiffness(fDimension*i+k,fDimension*nphi+j) += Q_SUPG * WJ;
                Stiffness(fDimension*nphi+j,fDimension*i+k) += Q * WJ;


                //PSPG stabilization
                double H = data.fDPhiX0(i,k) * shapeFj * tPSPG_ * fAlphaM;
                double G = data.fDPhiX0(i,k) * wSUPGj * tPSPG_ * AGDT;
                double Guu = 0.;
                for (int m = fDimension; m--; ) Guu += data.fDPhiX0(i,m)*duna_dx(m,k) * shapeFj * tPSPG_ * AGDT;

                Stiffness(fDimension*(nphi)+i,fDimension*j+k) += (H + G + Guu) * WJ;
            }

            double Q = 0.;
            for (int m = fDimension; m--; ) Q += data.fDPhiX0(i,m) * data.fDPhiX0(j,m) * tPSPG_ / fDensity;
            Stiffness(fDimension*nphi+i,fDimension*nphi+j) += Q * WJ;
        };
    };
}

void TransientNavierStokes::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    //Velocity
    VecDouble una_ = fAlphaF * data.fSol + (1. - fAlphaF) * data.fSolPrev;

    //Velocity Derivatives
    MatrixDouble duna_dx = fAlphaF * data.fDSolDx + (1. - fAlphaF) * data.fDSolDxPrev;
    
    //Mesh Velocity
    VecDouble umesh_(fDimension), umeshPrev_(fDimension), umeshna_(fDimension);
    umeshna_ = fAlphaF * umesh_ + (1. - fAlphaF) * umeshPrev_;

    //Pressure
    double p_ = data.fSol[fDimension];
    int nphi = data.fPhi.size();

    double AGDT = fAlphaF * fGamma * fTimeStep; 
    double VAGDT = fViscosity * AGDT; 
    double DAGDT = fDensity * AGDT; 
    double WJ = data.fWeight * data.fJacA0;
    double DAM = fDensity * fAlphaM;

    double divrU = 0.;
    for (int l=fDimension; l--; ) divrU += duna_dx(l,l);
    
    VecDouble forcingF(fDimension);
    forcingF.setZero();

    for (int i = nphi; i--; ){
        double shapeFi = data.fPhi[i];

        for (int k = fDimension; k--; ){

            //Mass + SUPG mass
            double m = shapeFi * data.fDSolDt[k] * fDensity;
            for (int l=fDimension; l--; ) m += (una_[l] - umeshna_[l]) * data.fDPhiX0(i,l) * data.fDSolDt[k] * tSUPG_ * fDensity;
                    
            //Viscosity
            double K = 0.;
            for (int l=fDimension; l--; ) K += data.fDPhiX0(i,l) * duna_dx(k,l) * fViscosity;
            for (int l=fDimension; l--; ) K += data.fDPhiX0(i,l) * duna_dx(l,k) * fViscosity;

            //LSIC
            double KLS = data.fDPhiX0(i,k) * divrU * tLSIC_ * fDensity;

            //Convection + SUPG
            double C = 0.;
            for (int l=fDimension; l--; ) C += duna_dx(k,l) * (una_[l] - umeshna_[l]) * shapeFi * fDensity;
            double conv = 0.;
            for (int l=fDimension; l--; ) conv += (una_[l] - umeshna_[l]) * data.fDPhiX0(i,l);
            for (int l=fDimension; l--; ) C += conv * (una_[l] - umeshna_[l]) * duna_dx(k,l) * tSUPG_ * fDensity;

            //Pressure + SUPG
            double P = - (data.fDPhiX0(i,k) * p_);
            for (int l=fDimension; l--; ) P += data.fDPhiX0(i,l) * (una_[l] - umeshna_[l]) * data.fDSolDx(fDimension,k) * tSUPG_;

            //External force
            double F = forcingF[k] * shapeFi * fDensity;
            
            Rhs[fDimension*i+k] += (-m -K - P - C - KLS + F) * WJ;
        }

        double Q = divrU * shapeFi;
        for (int l=fDimension; l--; ) Q += data.fDPhiX0(i,l) * data.fDSolDx(fDimension,l) * tPSPG_ / fDensity
                                  + data.fDPhiX0(i,l) * data.fDSolDt[l] * tPSPG_ 
                                  - data.fDPhiX0(i,l) * forcingF[l] * tPSPG_;
        for (int k=fDimension; k--; )
            for (int l=fDimension; l--; )
                Q += data.fDPhiX0(i,k) * (una_[l] - umeshna_[l]) * duna_dx(k,l) * tPSPG_;

        Rhs[fDimension*nphi+i] += -Q * WJ;
                            
    };
};

void TransientNavierStokes::ComputeError(IntPointData &data, VecDouble &errors){
    NavierStokes::ComputeError(data,errors);
}

int TransientNavierStokes::VariableIndex(const std::string &name) const{
    return NavierStokes::VariableIndex(name);
};

int TransientNavierStokes::NSolutionVariables(int var) const{
    return NavierStokes::NSolutionVariables(var);
};

void TransientNavierStokes::Solution(IntPointData &data, int var, VecDouble &Sol){
    NavierStokes::Solution(data,var,Sol);
};

void TransientNavierStokes::UpdateTimeDerivatives(CompMesh *cmesh){
    PanicButton();
};