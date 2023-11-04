#include "Stokes.h"

//------------------------------------------------------------------------------
//------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
//------------------------------------------------------------------------------
void Stokes::GetStabilizationParameter(int &index, IntPointData &data) {
    int DIM = fDimension;
    int nphi = data.fPhi.size();
    double tSUGN1_ = 0.;
    double tSUGN2_ = 0.;
    double tSUGN3_ = 0.;
    double hRGN_ = 0.;
    double r[DIM] = {};
    double s[DIM] = {};
    double hUGN_ = 0.;
    auto dphi_dx = data.fDPhiX0;
    tSUPG_ = 0.;  

    double u__[DIM] = {};
    double aux = 0.;
    double aux2 = 0.;

    for (int i = nphi; i--; ){
        double a1 = 0.;
        for (int j = DIM; j--; ){
            double ua = data.fSol(j);
            double uma = 0.;
            ua -= uma;
            u__[j] += ua * data.fPhi[i];
            a1 += ua*ua;
        }

        // a1 = std::sqrt(a1);
        for (int j = DIM; j--; ) r[j] += std::sqrt(a1) * dphi_dx(i,j);
    };


    for (int j = DIM; j--; ) aux += u__[j]*u__[j];
    double uNorm = std::sqrt(aux);

    // for (int i = fMesh->nElNodes; i--; ){
    //     for (int j = DIM; j--; ){
    //         r[j] += uNorm * dphi_dx[j][i];
    //     }
    // }

    for (int j = DIM; j--; ) aux2 += r[j]*r[j];
    double rNorm = std::sqrt(aux2);

    if(uNorm > 1.e-10){
        for (int j = DIM; j--; ) s[j] = u__[j] / uNorm;
    }else{
        for (int j = DIM; j--; ) s[j] = 1. / std::sqrt(2.);
    };

    if (rNorm >= 1.e-10){
        for (int j = DIM; j--; ) r[j] /= rNorm;
    }else{
        for (int j = DIM; j--; ) r[j] = 1. / std::sqrt(2.);
    };
    
    for (int i = nphi; i--; ){
        for (int j = DIM; j--; ){
            hRGN_ += r[j] * dphi_dx(i,j);
            hUGN_ += s[j] * dphi_dx(i,j);
        }
    };
    hRGN_ = std::fabs(hRGN_);
    hUGN_ = std::fabs(hUGN_);

    if (hRGN_ >= 1.e-10){
        hRGN_ = 2. / hRGN_;
    }else{
        hRGN_ = 2. / 1.e-10;
    };

    if (hUGN_ >= 1.e-10){
        hUGN_ = 2. / hUGN_;
    }else{
        hUGN_ = 2. / 1.e-10;
    };    

    if (uNorm >= 1.e-10){
        tSUGN1_ = hUGN_ / (2. * uNorm);
    }else{
        tSUGN1_ = hUGN_ / 2.e-10;
    };
              
    tSUGN2_ = 10000.;//dTime_ / 2.;

    tSUGN3_ = hRGN_ * hRGN_ / (4. * fViscosity / fDensity);
   
    if (std::fabs(tSUGN1_) <= 1.e-10) tSUGN1_ = 1.e-10;
    if (std::fabs(tSUGN3_) <= 1.e-10) tSUGN3_ = 1.e-10;

    //if(model == false) std::cout << "SUPG " << tSUGN1_ << " " << tSUGN3_ << std::endl;
    //Computing tSUPG parameter
    tSUPG_ = 1. / std::sqrt(1. / (tSUGN1_ * tSUGN1_) + 
                       1. / (tSUGN2_ * tSUGN2_) + 
                       1. / (tSUGN3_ * tSUGN3_));
    //tSUPG_ = 0.;


    tPSPG_ = 1.e-2*tSUPG_;
    tLSIC_ = tSUPG_ * uNorm * uNorm;


    return;
};

void Stokes::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    data.fNeedsDSol = true;
    data.fDSolDx.resize(fDimension, fNState);
    data.fNeedsSol = true;
    data.fSol.resize(fNState);
    int nphi = data.fPhi.size();

    int DIM = fDimension;
    auto dphi_dx = data.fDPhiX0;
    GetStabilizationParameter(index, data);

    // Trust me, it improves performance!
    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];

    for (int i = nphi; i-- ; ){        
        double shapeFi = data.fPhi[i];
        for (int j = nphi; j-- ; ){
            
            double shapeFj = data.fPhi[j];

            for (int k = DIM; k--;  ){
                for (int l = DIM; l--; ){

                    //Diffusion matrix
                    double K = dphi_dx(i,l) * dphi_dx(j,k) * fViscosity;
                    if (k==l) for (int m = DIM; m--; ) K += dphi_dx(i,m) * dphi_dx(j,m)* fViscosity;

                    Stiffness((DIM+1)*i+k,(DIM+1)*j+l) += K * WJ;
                }
                //Gradient operator
                //Divergent operator
                double Q = dphi_dx(i,k) * shapeFj;

                Stiffness((DIM+1)*i+k,(DIM+1)*j+DIM) += -Q * WJ;
                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += Q * WJ;
            }
            //PSPG stabilization
            double Q = 0.;
            for (int m = DIM; m--; ) Q += dphi_dx(i,m) * dphi_dx(j,m) * tPSPG_ / fDensity;
            Stiffness((DIM+1)*j+DIM,(DIM+1)*i+DIM) += Q * WJ;
        };
    };
    
}

void Stokes::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){
    int nphi = data.fPhi.size();
    auto force = fForceFunction;
    int DIM = fDimension;

    VecDouble forcingF(DIM);
    forcingF.setZero();
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    double WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    auto dphi_dx = data.fDPhiX0;

    double divrU = 0.;
    for (int l=DIM; l--; ) divrU += data.fDSolDx(l,l);

    for (int i = nphi; i--; ){
        // std::cout << "Sol = " << this->Mesh()->NodeVec()[this->getConnectivity()[i]]->GetSolution(0) << std::endl;
        double shapeFi = data.fPhi[i];

        for (int k = DIM; k--; ){
            //Viscosity
            double K = 0.;
            for (int l=DIM; l--; ) K += dphi_dx(i,l) * data.fDSolDx(k,l) * fViscosity;
            for (int l=DIM; l--; ) K += dphi_dx(i,l) * data.fDSolDx(l,k) * fViscosity;

            //Pressure + SUPG
            double P = - (dphi_dx(i,k) * data.fSol[DIM]);

            //External force
            double F = (forcingF[k]) * shapeFi;
            
            Rhs[(DIM+1)*i+k] += (-K - P + F) * WJ;
        }

        double Q = divrU * shapeFi;
        for (int l=DIM; l--; ) Q += dphi_dx(i,l) * data.fDSolDx(DIM,l) * tPSPG_ / fDensity
                                  + dphi_dx(i,l) * (forcingF[l]/fDensity) * tPSPG_;

        Rhs[(DIM+1)*i+DIM] += -Q * WJ;
                            
    };

};

void Stokes::ComputeError(IntPointData &data, VecDouble &errors){
    std::cout << "Not implemented yet\n";
    
    // PanicButton();
}
