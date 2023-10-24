#include "ElStokes.h"

//------------------------------------------------------------------------------
//------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
//------------------------------------------------------------------------------
void ElStokes::GetStabilizationParameter(int &index, double &tSUPG_, double &tPSPG_, double &tLSIC_, MatrixDouble &dphi_dx) {
    int DIM = Mesh()->Dimension();

    double tSUGN1_ = 0.;
    double tSUGN2_ = 0.;
    double tSUGN3_ = 0.;
    double hRGN_ = 0.;
    double r[DIM] = {};
    double s[DIM] = {};
    double hUGN_ = 0.;

    tSUPG_ = 0.;  

    double u__[DIM] = {};
    double aux = 0.;
    double aux2 = 0.;

    double &alpha_f = Mesh()->getProblemParameters().getAlphaF();
    double &visc_ = Mesh()->getProblemParameters().GetViscosity();
    double &dens_ = Mesh()->getProblemParameters().GetDensity();
    double &dTime_ = Mesh()->getProblemParameters().GetTimeStep();

    for (int i = Mesh()->NElNodes(); i--; ){
        double a1 = 0.;
        for (int j = DIM; j--; ){
            double ua =Mesh()->NodeVec()[getConnectivity()[i]] -> GetSolution(j);
            double uma = 0.;
            ua -= uma;
            u__[j] += ua * Mesh()->getNumericalIntegration()-> phi_(i,index);
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
    
    for (int i = Mesh()->NElNodes(); i--; ){
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
              
    tSUGN2_ = dTime_ / 2.;

    tSUGN3_ = hRGN_ * hRGN_ / (4. * visc_ / dens_);
   
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



void ElStokes::ComputeStiffness(int &index, MatrixDouble &Stiffness){

    double &visc_ = Mesh()->getProblemParameters().GetViscosity();
    double &dens_ = Mesh()->getProblemParameters().GetDensity();
    int DIM = Mesh()->Dimension();
    auto dphi_dx = fIntegData.fDPhiX0;
    GetStabilizationParameter(index, tSUPG_, tPSPG_, tLSIC_, dphi_dx);

    // Trust me, it improves performance!
    double WJ = fIntegData.fWeight * fIntegData.fJacA0 * getIntegPointWeightFunction(index);

    for (int i = Mesh()->NElNodes(); i-- ; ){        
        double shapeFi = Mesh()->getNumericalIntegration()-> phi_(i,index);
        for (int j = Mesh()->NElNodes(); j-- ; ){
            
            double shapeFj = Mesh()->getNumericalIntegration()-> phi_(j,index);

            for (int k = DIM; k--;  ){
                for (int l = DIM; l--; ){

                    //Diffusion matrix
                    double K = dphi_dx(i,l) * dphi_dx(j,k) * visc_;
                    if (k==l) for (int m = DIM; m--; ) K += dphi_dx(i,m) * dphi_dx(j,m)* visc_;

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
            for (int m = DIM; m--; ) Q += dphi_dx(i,m) * dphi_dx(j,m) * tPSPG_ / dens_;
            Stiffness((DIM+1)*j+DIM,(DIM+1)*i+DIM) += Q * WJ;
        };
    };
    
}

void ElStokes::ComputeResidual(int &index, VecDouble &Rhs){

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
    auto dphi_dx = fIntegData.fDPhiX0;

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

            //Pressure + SUPG
            double P = - (dphi_dx(i,k) * u_[DIM]);

            //External force
            double F = (fieldForce[k]*dens_ + forcingF[k]) * shapeFi;
            
            Rhs[(DIM+1)*i+k] += (-K - P + F) * WJ;
        }

        double Q = divrU * shapeFi;
        for (int l=DIM; l--; ) Q += dphi_dx(i,l) * du_dx(DIM,l) * tPSPG_ / dens_
                                  + dphi_dx(i,l) * (fieldForce[l] + forcingF[l]/dens_) * tPSPG_;

        Rhs[(DIM+1)*i+DIM] += -Q * WJ;
                            
    };

};

void ElStokes::ComputeError(VecDouble &errors){
    std::cout << "Not implemented yet\n";
    
    // PanicButton();
}

void ElStokes::ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs){
    int DIM = Mesh()->Dimension();
    // std::cout << "Stiffness antes = \n" << Stiffness<< std::endl;
    for (int i = Mesh()->NElNodes(); i--; ){
        int nstate = Mesh()->NodeVec()[getConnectivity()[i]]->GetNStateVariables();
        for (int istate = 0; istate < nstate-1; istate++){
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

        if (Mesh()->NodeVec()[getConnectivity()[i]]->getCoordinateValue(0) < 0.0001 &&
            Mesh()->NodeVec()[getConnectivity()[i]]->getCoordinateValue(1) < 0.0001){
            for (int j = Mesh()->NElNodes()*nstate; j--; ){
                Stiffness(nstate*i+DIM,j) = 0.;
                Stiffness(j,nstate*i+DIM) = 0.;
            };
            Stiffness(nstate*i+DIM,nstate*i+DIM) = 1.;
            Rhs[nstate*i+DIM] = 0.;
        }
    }
    // std::cout << "Stiffness depois = \n" << Stiffness<< std::endl;

}


