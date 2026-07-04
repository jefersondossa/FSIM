#include "TransientNavierStokes.h"

TransientNavierStokes::TransientNavierStokes(int matid, int dim, REAL density, REAL viscosity, REAL dTime, REAL specRadius) : NavierStokes(matid,dim,density,viscosity) {
    fTimeStep = dTime;
    fIntegScheme = EGeneralizedAlpha;
    fSpectralRadius = specRadius;
    fAlphaF = 1. / (1. + specRadius);
    fAlphaM = 0.5 * (3. - specRadius)/(1. + specRadius);
    fGamma = 0.5 + fAlphaM - fAlphaF;
    this->fType = WeakFormType::kTransientNavierStokes;
};

//Stabilization parameters computation
void TransientNavierStokes::GetStabilizationParameter(int &index, IntPointData &data) {
    int DIM = fDimension;
    int nphi = data.fPhi.size();
    REAL tSUGN1_ = 0.;
    REAL tSUGN2_ = 0.;
    REAL tSUGN3_ = 0.;
    REAL hRGN_ = 0.;
    REAL r[DIM] = {};
    REAL s[DIM] = {};
    REAL hUGN_ = 0.;
    auto dphi_dx = data.fDPhiX0;
    tSUPG_ = 0.;  

    REAL u__[DIM] = {};
    REAL aux = 0.;
    REAL aux2 = 0.;

    for (int i = nphi; i--; ){
        REAL a1 = 0.;
        for (int j = DIM; j--; ){
            REAL ua = data.fSol(j);
            REAL uma = 0.;
            ua -= uma;
            u__[j] += ua * data.fPhi[i];
            a1 += ua*ua;
        }

        // a1 = std::sqrt(a1);
        for (int j = DIM; j--; ) r[j] += std::sqrt(a1) * dphi_dx(j,i);
    };


    for (int j = DIM; j--; ) aux += u__[j]*u__[j];
    REAL uNorm = std::sqrt(aux);

    // for (int i = fMesh->nElNodes; i--; ){
    //     for (int j = DIM; j--; ){
    //         r[j] += uNorm * dphi_dx[j][i];
    //     }
    // }

    for (int j = DIM; j--; ) aux2 += r[j]*r[j];
    REAL rNorm = std::sqrt(aux2);

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
            hRGN_ += r[j] * dphi_dx(j,i);
            hUGN_ += s[j] * dphi_dx(j,i);
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
              
    tSUGN2_ = fTimeStep / 2.;

    tSUGN3_ = hRGN_ * hRGN_ / (4. * fViscosity / fDensity);
   
    if (std::fabs(tSUGN1_) <= 1.e-10) tSUGN1_ = 1.e-10;
    if (std::fabs(tSUGN3_) <= 1.e-10) tSUGN3_ = 1.e-10;

    //if(model == false) std::cout << "SUPG " << tSUGN1_ << " " << tSUGN3_ << std::endl;
    //Computing tSUPG parameter
    tSUPG_ = 1. / std::sqrt(1. / (tSUGN1_ * tSUGN1_) + 
                       1. / (tSUGN2_ * tSUGN2_) + 
                       1. / (tSUGN3_ * tSUGN3_));
    //tSUPG_ = 0.;


    tPSPG_ = tSUPG_;
    tLSIC_ = tSUPG_ * uNorm * uNorm;
 

    return;
};


void TransientNavierStokes::ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness){
    
    // if (fTimeStep < 1){
    //     NavierStokes::ComputeStiffness(index,data,Stiffness);
    //     return;
    // }

    GetStabilizationParameter(index, data);
    
    // NavierStokes::ComputeStiffness(index,data,Stiffness);
    // Stiffness *= fAlphaF;
    //Velocity
    VecDouble una_ = fAlphaF * data.fSol + (1. - fAlphaF) * data.fSolPrev;

    //Velocity Derivatives
    MatrixDouble duna_dx = fAlphaF * data.fDSolDx + (1. - fAlphaF) * data.fDSolDxPrev;
    
    //Mesh Velocity
    VecDouble umesh_(fDimension), umeshPrev_(fDimension), umeshna_(fDimension);
    umeshna_ = fAlphaF * umesh_ + (1. - fAlphaF) * umeshPrev_;

    // // Trust me, it improves performance!
    REAL AGDT = fAlphaF; 
    // REAL VAGDT = fViscosity * AGDT; 
    // REAL DAGDT = fDensity * AGDT; 
    // REAL WJ = data.fWeight * data.fJacA0;
    REAL DAM = fDensity * fAlphaM / (fGamma * fTimeStep);

    // tSUPG_ = 0.;
    // tLSIC_ = 0.;

    // int nphi = data.fPhi.size();
    // for (int i = nphi; i-- ; ){
        
    //     REAL shapeFi = data.fPhi[i];
    //     REAL wSUPGi = 0.;
    //     for (int m=fDimension; m--; ) wSUPGi += (una_[m] - umeshna_[m]) * data.fDPhiX0(m,i);

    //     for (int j = nphi; j-- ; ){
            
    //         REAL shapeFj = data.fPhi[j];
    //         REAL shapeFij = shapeFi * shapeFj;
    //         REAL wSUPGj = 0.;
    //         for (int m=fDimension; m--; ) wSUPGj += (una_[m] - umeshna_[m]) * data.fDPhiX0(m,j)*0.;
            
    //         //Mass matrix (use for both directions)
    //         REAL M = (shapeFij + wSUPGi * shapeFj * tSUPG_) * DAM*0.;

    //         //Convection matrix
    //         REAL C = (wSUPGj * shapeFi + wSUPGi * wSUPGj * tSUPG_) * DAGDT;         

    //         REAL aux1 = tSUPG_ * wSUPGi * shapeFj;
    //         REAL aux2 = tSUPG_ * shapeFj;

    //         for (int k = fDimension; k--;  ){

    //             Stiffness(fDimension*i+k,fDimension*j+k) += (M + C) * WJ;
                    
    //             REAL conv = 0.;
    //             for (int m = fDimension; m--; ) conv += una_[m]*duna_dx(k,m)*0.;

    //             for (int l = fDimension; l--; ){

    //                 //Diffusion matrix
    //                 REAL K = data.fDPhiX0(l,i) * data.fDPhiX0(k,j) * VAGDT;
    //                 if (k==l) for (int m = fDimension; m--; ) K += data.fDPhiX0(m,i) * data.fDPhiX0(m,j) * VAGDT;

    //                 //Convection derivatives
    //                 REAL Cuu = (shapeFij * duna_dx(k,l) + 
    //                               aux1 * duna_dx(k,l) +
    //                               aux2 * conv * data.fDPhiX0(l,i)) * DAGDT*0.;

    //                 //LSIC
    //                 REAL KLS = data.fDPhiX0(k,i) * data.fDPhiX0(l,j) * tLSIC_ * DAGDT*0.;

    //                 Stiffness(fDimension*i+k,fDimension*j+l) += (K + KLS + Cuu) * WJ;
    //             }

    //             //SINAL DA PARCELA QUE MULTIPLICA O TSUPG ESTA COM SINAL TROCADO NA FORMULAÇAO DO TEZDUYAR
    //             //Gradient operator
    //             REAL Q_SUPG = - data.fDPhiX0(k,i) * shapeFj + wSUPGi * data.fDPhiX0(k,j) * tSUPG_;
    //             //Divergent operator
    //             REAL Q = data.fDPhiX0(k,i) * shapeFj * AGDT;

    //             Stiffness(fDimension*i+k,fDimension*nphi+j) += Q_SUPG * WJ;
    //             Stiffness(fDimension*nphi+j,fDimension*i+k) += Q * WJ;


    //             //PSPG stabilization
    //             REAL H = data.fDPhiX0(k,i) * shapeFj * tPSPG_ * fAlphaM*0.;
    //             REAL G = data.fDPhiX0(k,i) * wSUPGj * tPSPG_ * AGDT*0.;
    //             REAL Guu = 0.;
    //             for (int m = fDimension; m--; ) Guu += data.fDPhiX0(m,i)*duna_dx(m,k) * shapeFj * tPSPG_ * AGDT*0.;

    //             Stiffness(fDimension*(nphi)+i,fDimension*j+k) += (H + G + Guu) * WJ;
    //         }

    //         REAL Q = 0.;
    //         for (int m = fDimension; m--; ) Q += data.fDPhiX0(m,i) * data.fDPhiX0(m,j) * tPSPG_ / fDensity;
    //         Stiffness(fDimension*nphi+i,fDimension*nphi+j) += Q * WJ;
    //     };
    // };

    int DIM = fDimension;
    auto dphi_dx = data.fDPhiX0;
    // this->GetStabilizationParameter(index, data);
    REAL WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];
    // this->tSUPG_ = 0.;
    // this->tLSIC_ = 0.;
    int nphi = data.fPhi.size();

    for (int i = nphi; i-- ; ){        
        REAL shapeFi = data.fPhi[i];
        REAL wSUPGi = 0.;
        for (int m=DIM; m--; ) wSUPGi += una_[m] * dphi_dx(m,i);

        for (int j = nphi; j-- ; ){
            
            REAL shapeFj = data.fPhi[j];
            REAL shapeFij = shapeFi * shapeFj;
            REAL wSUPGj = 0.;
            for (int m=DIM; m--; ) wSUPGj += una_[m] * dphi_dx(m,j);

            //Mass matrix (use for both directions)
            REAL M = (shapeFij + wSUPGi * shapeFj * tSUPG_) * DAM;

            //Convection matrix
            REAL C = (wSUPGj * shapeFi + wSUPGi * wSUPGj * this->tSUPG_) * fDensity * fAlphaF;

            REAL aux1 = this->tSUPG_ * wSUPGi * shapeFj;
            REAL aux2 = this->tSUPG_ * shapeFj;

            for (int k = DIM; k--;  ){

                Stiffness((DIM+1)*i+k,(DIM+1)*j+k) += (M+C) * WJ;
                
                REAL conv = 0.;
                for (int m = DIM; m--; ) conv += una_[m]*duna_dx(k,m);

                for (int l = DIM; l--; ){

                    //Diffusion matrix
                    REAL K = dphi_dx(l,i) * dphi_dx(k,j) * fViscosity;
                    if (k==l) for (int m = DIM; m--; ) K += dphi_dx(m,i) * dphi_dx(m,j)* fViscosity;

                    //Convection derivatives
                    REAL Cuu = (shapeFij * duna_dx(k,l) + 
                                  aux1 * duna_dx(k,l) +
                                  aux2 * conv * data.fDPhiX0(l,i)) * fDensity;

                    //LSIC
                    REAL KLS = data.fDPhiX0(k,i) * data.fDPhiX0(l,j) * this->tLSIC_ * fDensity;

                    Stiffness((DIM+1)*i+k,(DIM+1)*j+l) += (KLS + Cuu + K) * WJ * fAlphaF;
                }
                //Gradient operator
                //Divergent operator
                REAL Q = dphi_dx(k,i) * shapeFj;

                Stiffness((DIM+1)*i+k,(DIM+1)*j+DIM) += -Q * WJ;
                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += Q * WJ * fAlphaF;

                //Gradient operator
                REAL Q_SUPG = wSUPGi * data.fDPhiX0(k,j) * this->tSUPG_;

                Stiffness((DIM+1)*i+k,(DIM+1)*j+DIM) += Q_SUPG * WJ;

                //PSPG stabilization
                REAL G = data.fDPhiX0(k,i) * wSUPGj * this->tPSPG_ * fAlphaF;
                REAL H = data.fDPhiX0(k,i) * shapeFj * tPSPG_ * DAM;

                REAL Guu = 0.;
                for (int m = DIM; m--; ) Guu += data.fDPhiX0(m,i) * duna_dx(m,k) * shapeFj * this->tPSPG_ * fAlphaF;

                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += (G + Guu + H) * WJ;
            }
            //PSPG stabilization
            REAL Q = 0.;
            for (int m = DIM; m--; ) Q += dphi_dx(m,i) * dphi_dx(m,j) * tPSPG_ / fDensity;
            Stiffness((DIM+1)*j+DIM,(DIM+1)*i+DIM) += Q * WJ;
        };
    };

}

void TransientNavierStokes::ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs){

    // if (fTimeStep < 1 ){
    //     NavierStokes::ComputeResidual(index,data,Rhs);
    //     return;
    // }
    
    auto force = fForceFunction;
    int DIM = fDimension;

    int nphi = data.fPhi.size();
    VecDouble forcingF(DIM);
    forcingF.setZero();
    VecDouble x_ = data.fX;
    if (force) force(x_,forcingF);

    //Velocity
    VecDouble una_ = fAlphaF * data.fSol + (1. - fAlphaF) * data.fSolPrev;

    //Velocity Derivatives
    MatrixDouble duna_dx = fAlphaF * data.fDSolDx + (1. - fAlphaF) * data.fDSolDxPrev;
    
    //Mesh Velocity
    VecDouble umesh_(fDimension), umeshPrev_(fDimension), umeshna_(fDimension);
    umeshna_ = fAlphaF * umesh_ + (1. - fAlphaF) * umeshPrev_;

    //Pressure
    REAL p_ = data.fSol[fDimension];

    REAL WJ = data.fWeight * data.fJacA0 * data.fWeightFunction[index];

    REAL divrU = 0.;
    for (int l=DIM; l--; ) divrU += duna_dx(l,l);

    // for (int i = nphi; i--; ){
    //     REAL shapeFi = data.fPhi[i];

    //     for (int k = DIM; k--; ){
    //         //Mass + SUPG mass
    //         REAL m = shapeFi * data.fDSolDt[k] * fDensity;
    //         for (int l=fDimension; l--; ) m += (una_[l] - umeshna_[l]) * data.fDPhiX0(l,i) * data.fDSolDt[k] * tSUPG_ * fDensity;

    //         //Viscosity
    //         REAL K = 0.;
    //         for (int l=DIM; l--; ) K += data.fDPhiX0(l,i) * duna_dx(k,l) * fViscosity;
    //         for (int l=DIM; l--; ) K += data.fDPhiX0(l,i) * duna_dx(l,k) * fViscosity;

    //         //LSIC
    //         REAL KLS = data.fDPhiX0(k,i) * divrU * this->tLSIC_ * fDensity;

    //         //Convection + SUPG
    //         REAL C = 0.;
    //         for (int l=DIM; l--; ) C += duna_dx(k,l) * una_[l] * shapeFi * fDensity;
    //         REAL conv = 0.;
    //         for (int l=DIM; l--; ) conv += una_[l] * data.fDPhiX0(l,i);
    //         for (int l=DIM; l--; ) C += conv * una_[l] * duna_dx(k,l) * this->tSUPG_ * fDensity;

    //         //Pressure + SUPG
    //         REAL P = - (data.fDPhiX0(k,i) * una_[DIM]);
    //         for (int l=DIM; l--; ) P += data.fDPhiX0(l,i) * una_[l] * duna_dx(DIM,k) * this->tSUPG_;

    //         //External force
    //         REAL F = (forcingF[k]) * shapeFi;
            
    //         Rhs[(DIM+1)*i+k] += (-m -K - P - C - KLS + F) * WJ;
    //     }

    //     REAL Q = divrU * shapeFi;
    //     for (int l=DIM; l--; ) Q += data.fDPhiX0(l,i) * duna_dx(DIM,l) * this->tPSPG_ / fDensity
    //                               + data.fDPhiX0(l,i) * data.fDSolDt[l] * tPSPG_
    //                               - data.fDPhiX0(l,i) * (forcingF[l]/fDensity) * this->tPSPG_;
    //     for (int k=DIM; k--; )
    //         for (int l=DIM; l--; )
    //             Q += data.fDPhiX0(k,i) * una_[l] * duna_dx(k,l) * this->tPSPG_;

    //     Rhs[(DIM+1)*i+DIM] += -Q * WJ;
                            
    // };




    for (int i = nphi; i--; ){
        REAL shapeFi = data.fPhi[i];

        for (int k = fDimension; k--; ){

            //Mass + SUPG mass
            REAL m = shapeFi * data.fDSolDt[k] * fDensity;
            for (int l=fDimension; l--; ) m += (una_[l] - umeshna_[l]) * data.fDPhiX0(l,i) * data.fDSolDt[k] * tSUPG_ * fDensity;
                    
            //Viscosity
            REAL K = 0.;
            for (int l=fDimension; l--; ) K += data.fDPhiX0(l,i) * duna_dx(k,l) * fViscosity;
            for (int l=fDimension; l--; ) K += data.fDPhiX0(l,i) * duna_dx(l,k) * fViscosity;

            //LSIC
            REAL KLS = data.fDPhiX0(k,i) * divrU * tLSIC_ * fDensity;

            //Convection + SUPG
            REAL C = 0.;
            for (int l=fDimension; l--; ) C += duna_dx(k,l) * (una_[l] - umeshna_[l]) * shapeFi * fDensity;
            REAL conv = 0.;
            for (int l=fDimension; l--; ) conv += (una_[l] - umeshna_[l]) * data.fDPhiX0(l,i);
            for (int l=fDimension; l--; ) C += conv * (una_[l] - umeshna_[l]) * duna_dx(k,l) * tSUPG_ * fDensity;

            //Pressure + SUPG
            REAL P = - (data.fDPhiX0(k,i) * p_);
            for (int l=fDimension; l--; ) P += data.fDPhiX0(l,i) * (una_[l] - umeshna_[l]) * data.fDSolDx(fDimension,k) * tSUPG_;

            //External force
            REAL F = forcingF[k] * shapeFi * fDensity;
            
            Rhs[fDimension*i+k] += (-m -K - P - C - KLS + F) * WJ;
        }

        REAL Q = divrU * shapeFi;
        for (int l=fDimension; l--; ) Q += data.fDPhiX0(l,i) * data.fDSolDx(fDimension,l) * tPSPG_ / fDensity
                                  + data.fDPhiX0(l,i) * data.fDSolDt[l] * tPSPG_ 
                                  - data.fDPhiX0(l,i) * forcingF[l] * tPSPG_;
        for (int k=fDimension; k--; )
            for (int l=fDimension; l--; )
                Q += data.fDPhiX0(k,i) * (una_[l] - umeshna_[l]) * duna_dx(k,l) * tPSPG_;

        Rhs[fDimension*nphi+i] += -Q * WJ;
                            
    };
};

void TransientNavierStokes::ComputeError(IntPointData &data, VecDouble &errors){
    NavierStokes::ComputeError(data,errors);
}

int TransientNavierStokes::VariableIndex(const std::string &name) const{
    int var = NavierStokes::VariableIndex(name);

    if (var != -1){
        return var;
    } else {
        if(!strcmp("Acceleration",name.c_str()))       return 11;
    } 

    return -1;
};

int TransientNavierStokes::NSolutionVariables(int var) const{
    int nsol = NavierStokes::NSolutionVariables(var);
    if (nsol != -1){
        return nsol;
    } else {
        if (var == 11) return 3;
    }
    return -1;
};

void TransientNavierStokes::Solution(IntPointData &data, int var, VecDouble &Sol){
    NavierStokes::Solution(data,var,Sol);
    //Acceleration
    if (var == 11){
        Sol[0] = data.fDSolDt[0];
        Sol[1] = data.fDSolDt[1];
        Sol[2] = 0.;
        return;
    };
};

void TransientNavierStokes::UpdateTimeDerivatives(CompMesh *cmesh){
    
    for (int inode = 0; inode < cmesh->NConnects(); inode++){
        for (int j = fDimension; j--; ){
            REAL u_ = cmesh->ConnectVec()[inode] -> GetSolution(j);
            REAL uprev = cmesh->ConnectVec()[inode] -> GetPreviousSolution(j);
            // cmesh->NodeVec()[inode] -> SetPreviousSolution(j,u_);

            REAL accel = cmesh->ConnectVec()[inode] -> GetDSolutionDTime(j);
            // cmesh->NodeVec()[inode] -> setPreviousAccelerationComponent(j,accel);
             
            REAL accelUpdated = (u_- uprev)/(fGamma * fTimeStep) + accel * (fGamma - 1.) / fGamma;
            cmesh->ConnectVec()[inode] -> SetDSolutionDTime(j,accelUpdated);

            // int dof_i = (fDimension+1)*inode + j;
            // VecSetValues(u, 1, &dof_i, &accel,INSERT_VALUES);
        }
    }
};