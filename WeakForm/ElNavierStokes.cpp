#include "ElNavierStokes.h"

template<class tshape>
void ElNavierStokes<tshape>::ComputeStiffness(int &index, MatrixDouble &Stiffness){
    double &visc_ = this->Mesh()->getProblemParameters().GetViscosity();
    double &dens_ = this->Mesh()->getProblemParameters().GetDensity();
    int DIM = this->Mesh()->Dimension();

    this->GetStabilizationParameter(index, this->tSUPG_, this->tPSPG_, this->tLSIC_, this->fIntegData.fDPhiX0);
    double WJ = this->fIntegData.fWeight * this->fIntegData.fJacA0 * this->getIntegPointWeightFunction(index);
    this->tSUPG_ = 0.;
    this->tLSIC_ = 0.;

    VecDouble u_(DIM+1);
    this->interpolateSolution(index, u_);

    //Solution Derivatives
    MatrixDouble du_dx(DIM+1,DIM);
    this->interpolateSolDerivatives(du_dx);
    auto dphi_dx = this->fIntegData.fDPhiX0;

    for (int i = tshape::NElNodes; i-- ; ){        
        double shapeFi = this->fIntegData.fPhi[i];
        double wSUPGi = 0.;
        for (int m=DIM; m--; ) wSUPGi += u_[m] * dphi_dx(i,m);
 
        for (int j = tshape::NElNodes; j-- ; ){           

            double shapeFj = this->fIntegData.fPhi[j];
            double shapeFij = shapeFi * shapeFj;
            double wSUPGj = 0.;
            for (int m=DIM; m--; ) wSUPGj += u_[m] * dphi_dx(j,m);
            
            //Convection matrix
            double C = (wSUPGj * shapeFi + wSUPGi * wSUPGj * this->tSUPG_) * dens_ * 0.;

            double aux1 = this->tSUPG_ * wSUPGi * shapeFj;
            double aux2 = this->tSUPG_ * shapeFj;

            for (int k = DIM; k--;  ){
                Stiffness((DIM+1)*i+k,(DIM+1)*j+k) += C * WJ;
                double conv = 0.;
                for (int m = DIM; m--; ) conv += u_[m]*du_dx(k,m);

                for (int l = DIM; l--; ){

                    //Diffusion matrix
                    double K = this->fIntegData.fDPhiX0(i,l) * this->fIntegData.fDPhiX0(j,k) * visc_;
                    if (k==l) for (int m = DIM; m--; ) K += this->fIntegData.fDPhiX0(i,m) * this->fIntegData.fDPhiX0(j,m)* visc_;

                    //Convection derivatives
                    double Cuu = (shapeFij * du_dx(k,l) + 
                                  aux1 * du_dx(k,l) +
                                  aux2 * conv * this->fIntegData.fDPhiX0(i,l)) * dens_;

                    //LSIC
                    double KLS = this->fIntegData.fDPhiX0(i,k) * this->fIntegData.fDPhiX0(j,l) * this->tLSIC_ * dens_;

                    Stiffness((DIM+1)*i+k,(DIM+1)*j+l) += (K + KLS + Cuu) * WJ;
                }
                //Gradient operator
                double Q_SUPG = - this->fIntegData.fDPhiX0(i,k) * shapeFj + wSUPGi * this->fIntegData.fDPhiX0(j,k) * this->tSUPG_;
                //Divergent operator
                double Q = this->fIntegData.fDPhiX0(i,k) * shapeFj;

                Stiffness((DIM+1)*i+k,(DIM+1)*j+DIM) += Q_SUPG * WJ;
                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += Q * WJ;

                //PSPG stabilization
                double G = this->fIntegData.fDPhiX0(i,k) * wSUPGj * this->tPSPG_;
                double Guu = 0.;
                for (int m = DIM; m--; ) Guu += this->fIntegData.fDPhiX0(i,m) * du_dx(m,k) * shapeFj * this->tPSPG_;

                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += (G + Guu) * WJ;
            }
            //PSPG stabilization
            double Q = 0.;
            for (int m = DIM; m--; ) Q += this->fIntegData.fDPhiX0(i,m) * this->fIntegData.fDPhiX0(j,m) * this->tPSPG_ / dens_;
            Stiffness((DIM+1)*j+DIM,(DIM+1)*i+DIM) += Q * WJ;
        };
    };



    for (int i = tshape::NElNodes; i-- ; ){        
        double shapeFi = this->fIntegData.fPhi[i];
        for (int j = tshape::NElNodes; j-- ; ){
            
            double shapeFj = this->fIntegData.fPhi[j];

            for (int k = DIM; k--;  ){
                for (int l = DIM; l--; ){

                    //Diffusion matrix
                    double K = this->fIntegData.fDPhiX0(i,l) * this->fIntegData.fDPhiX0(j,k) * visc_;
                    if (k==l) for (int m = DIM; m--; ) K += this->fIntegData.fDPhiX0(i,m) * this->fIntegData.fDPhiX0(j,m)* visc_;

                    Stiffness((DIM+1)*i+k,(DIM+1)*j+l) += K * WJ;
                }
                //Gradient operator
                //Divergent operator
                double Q = this->fIntegData.fDPhiX0(i,k) * shapeFj;

                Stiffness((DIM+1)*i+k,(DIM+1)*j+DIM) += -Q * WJ;
                Stiffness((DIM+1)*j+DIM,(DIM+1)*i+k) += Q * WJ;
            }
            //PSPG stabilization
            double Q = 0.;
            for (int m = DIM; m--; ) Q += this->fIntegData.fDPhiX0(i,m) * this->fIntegData.fDPhiX0(j,m) * this->tPSPG_ / dens_;
            Stiffness((DIM+1)*j+DIM,(DIM+1)*i+DIM) += Q * WJ;
        };
    };

}

template<class tshape>
void ElNavierStokes<tshape>::ComputeResidual(int &index, VecDouble &Rhs){
    VecDouble fieldForce = this->Mesh()->getProblemParameters().GetFieldForce();
    auto force = this->Mesh()->getProblemParameters().getForcingFunction();
    int DIM = this->Mesh()->Dimension();
    double &visc_ = this->Mesh()->getProblemParameters().GetViscosity();
    double &dens_ = this->Mesh()->getProblemParameters().GetDensity();

    VecDouble forcingF(DIM);
    forcingF.setZero();
    VecDouble x_ = this->getIntegPointCoordinatesValue(index);
    if (force) force(x_,forcingF);

    //Solution Derivatives
    MatrixDouble du_dx(DIM+1,DIM);
    this->interpolateSolDerivatives(du_dx);

    VecDouble u_(DIM+1);
    this->interpolateSolution(index, u_);

    double WJ = this->fIntegData.fWeight * this->fIntegData.fJacA0 * this->getIntegPointWeightFunction(index);

    double divrU = 0.;
    for (int l=DIM; l--; ) divrU += du_dx(l,l);

    for (int i = tshape::NElNodes; i--; ){
        // std::cout << "Sol = " << this->Mesh()->NodeVec()[this->getConnectivity()[i]]->GetSolution(0) << std::endl;
        double shapeFi = this->fIntegData.fPhi[i];

        for (int k = DIM; k--; ){
            //Viscosity
            double K = 0.;
            for (int l=DIM; l--; ) K += this->fIntegData.fDPhiX0(i,l) * du_dx(k,l) * visc_;
            for (int l=DIM; l--; ) K += this->fIntegData.fDPhiX0(i,l) * du_dx(l,k) * visc_;

            //LSIC
            double KLS = this->fIntegData.fDPhiX0(i,k) * divrU * this->tLSIC_ * dens_;

            //Convection + SUPG
            double C = 0.;
            for (int l=DIM; l--; ) C += du_dx(k,l) * u_[l] * shapeFi * dens_;
            double conv = 0.;
            for (int l=DIM; l--; ) conv += u_[l] * this->fIntegData.fDPhiX0(i,l);
            for (int l=DIM; l--; ) C += conv * u_[l] * du_dx(k,l) * this->tSUPG_ * dens_;

            //Pressure + SUPG
            double P = - (this->fIntegData.fDPhiX0(i,k) * u_[DIM]);
            for (int l=DIM; l--; ) P += this->fIntegData.fDPhiX0(i,l) * u_[l] * du_dx(DIM,k) * this->tSUPG_;

            //External force
            double F = (fieldForce[k]*dens_ + forcingF[k]) * shapeFi;
            
            Rhs[(DIM+1)*i+k] += (-K - P - C - KLS + F) * WJ;
        }

        double Q = divrU * shapeFi;
        for (int l=DIM; l--; ) Q += this->fIntegData.fDPhiX0(i,l) * du_dx(DIM,l) * this->tPSPG_ / dens_
                                  + this->fIntegData.fDPhiX0(i,l) * (fieldForce[l] + forcingF[l]/dens_) * this->tPSPG_;
        for (int k=DIM; k--; )
            for (int l=DIM; l--; )
                Q += this->fIntegData.fDPhiX0(i,k) * u_[l] * du_dx(k,l) * this->tPSPG_;

        Rhs[(DIM+1)*i+DIM] += -Q * WJ;
                            
    };
    
};

template<class tshape>
void ElNavierStokes<tshape>::ComputeError(VecDouble &errors){
    std::cout << "Not implemented yet.\n";
    // PanicButton();
}


#include "ShapeHexahedron.h"
#include "ShapeOneD.h"
#include "ShapeQuadrilateralLin.h"
#include "ShapePoint.h"
#include "ShapeTetrahedronLin.h"
#include "ShapeTetrahedronQua.h"
#include "ShapeTetrahedronCub.h"
#include "ShapeTriangleLin.h"
#include "ShapeTriangleQua.h"
#include "ShapeTriangleCub.h"

template class ElNavierStokes<ShapePoint>;
template class ElNavierStokes<ShapeOneD>;
template class ElNavierStokes<ShapeTriangleLin>;
template class ElNavierStokes<ShapeTriangleQua>;
template class ElNavierStokes<ShapeTriangleCub>;
template class ElNavierStokes<ShapeQuadrilateralLin>;
template class ElNavierStokes<ShapeTetrahedronLin>;
template class ElNavierStokes<ShapeTetrahedronQua>;
template class ElNavierStokes<ShapeTetrahedronCub>;
template class ElNavierStokes<ShapeHexahedron>;