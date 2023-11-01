#include "ElCouplingLocal.h"

template<class tshape>
void ElCouplingLocal<tshape>::ComputeStiffness(int &index, std::vector<MatrixDouble> &Stiffness){

    //Stiffness[0] is always the L2-H1 contribution;
    //Stiffness[1] is always the diagonal contribution to the global stiffness matrix
    //The remaining positions corresponds to the other contributions, 
    //which will be added to stiffness[0] in the global matrix and depends on the problem.


    int DIM = tshape::Dimension;
    int DEG = this->Mesh()->GetDefaultOrder();

    double k1 = fMeshVector[1]->getProblemParameters().getArlequinK1();
    double k2 = fMeshVector[1]->getProblemParameters().getArlequinK2();
        
    double WJ = this->fIntegData.fWeight * this->fIntegData.fJacA0;

    if (this->Mesh()->getProblemParameters().ProbType() == ProblemType::EPoisson){
        for (int i = 0; i < this->Mesh()->NElNodes(); i++){
            for (int j = 0; j < this->Mesh()->NElNodes(); j++){
                double l2 = this->Mesh()->getNumericalIntegration()-> phi_(i,index) * this->Mesh()->getNumericalIntegration()-> phi_(j,index) * WJ * k1;
                // L2 COUPLING OPERATOR
                Stiffness[0](i,j) -= l2;
                for (int l = 0; l < DIM; l++){
                    //H1 COUPLING OPERATOR
                    double K = this->fIntegData.fDPhiX0(i,l) * this->fIntegData.fDPhiX0(j,l);
                    Stiffness[0](i,j) -= K * WJ * k2;
                };
            };
        };
    } else {
        for (int i = 0; i < this->Mesh()->NElNodes(); i++){
            for (int j = 0; j < this->Mesh()->NElNodes(); j++){
                double l2 = this->Mesh()->getNumericalIntegration()-> phi_(i,index) * this->Mesh()->getNumericalIntegration()-> phi_(j,index) * WJ * k1;
                for (int k = 0; k < DIM; k++){
                    // L2 COUPLING OPERATOR
                    Stiffness[0](DIM*i+k,DIM*j+k) -= l2;
                    for (int l = 0; l < DIM; l++){
                        //H1 COUPLING OPERATOR
                        double K = this->fIntegData.fDPhiX0(i,l) * this->fIntegData.fDPhiX0(j,k);
                        if (k==l) for (int m = DIM; m--; ) K += this->fIntegData.fDPhiX0(i,m) * this->fIntegData.fDPhiX0(j,m);

                        Stiffness[0](DIM*i+k,DIM*j+l) -= K * WJ * k2;
                    };//l
                };//k
            };//j
        };//i
    };//if Poisson

    if (this->Mesh()->getProblemParameters().ArlequinStab() != ENoStab){
        ArlequinStabStiffness(index, this->fIntegData.fDPhiX0, this->fIntegData.fWeight, this->fIntegData.fJacA0, Stiffness);
    }

};

template<class tshape>
void ElCouplingLocal<tshape>::ComputeResidual(int &index, std::vector<VecDouble> &Rhs){

    int DIM = tshape::Dimension;
    int DEG = this->Mesh()->GetDefaultOrder();

    auto ellocal = fMeshVector[1]->ElementVec()[fLocalIndex];

    //Velocity
    VecDouble u_(this->Mesh()->NState()+1);
    ellocal->interpolateSolution(index, u_);

    //Velocity Derivatives
    MatrixDouble du_dx(this->Mesh()->NState()+1,DIM);
    ellocal->interpolateSolDerivatives(du_dx);

    //The lagrange multiplier is the solution of meshvector[2]
    //Lagrange Multiplier
    VecDouble lagM_(this->Mesh()->NState()+1);
    this->interpolateSolution(index, lagM_);

    //Lagrange Multiplier Derivatives
    MatrixDouble dL_dx(this->Mesh()->NState()+1,DIM);
    this->interpolateSolDerivatives(dL_dx);

    double k1 = fMeshVector[1]->getProblemParameters().getArlequinK1();
    double k2 = fMeshVector[1]->getProblemParameters().getArlequinK2();
        
    double WJ = this->fIntegData.fWeight * this->fIntegData.fJacA0;

    if (this->Mesh()->getProblemParameters().ProbType() == ProblemType::EPoisson){
        for (int i = 0; i < this->Mesh()->NElNodes(); i++){
            //Solution residual
            double L2u = u_[0] * this->Mesh()->getNumericalIntegration()-> phi_(i,index) * k1;
            double H1u = 0.;
            for (int l=DIM; l--; ) H1u += this->fIntegData.fDPhiX0(i,l) * du_dx(0,l) * k2;
            Rhs[0][this->Mesh()->NLocDOF()+i] += (L2u + H1u) * WJ;
            
            // Lagrange multipliers residual
            double L2 = lagM_[0] * this->Mesh()->getNumericalIntegration()-> phi_(i,index) * k1;
            double H1 = 0.;
            for (int l=DIM; l--; ) H1 += this->fIntegData.fDPhiX0(i,l) * dL_dx(0,l) * k2;
            Rhs[0][i] += (L2 + H1) * WJ;
        };
    } else {
        for (int i = 0; i < this->Mesh()->NElNodes(); i++){
            for (int k = 0; k < DIM; k++){
                //Solution residual
                double L2u = u_[k] * this->Mesh()->getNumericalIntegration()-> phi_(i,index) * k1;
                double H1u = 0.;
                for (int l=DIM; l--; ) H1u += this->fIntegData.fDPhiX0(i,l) * du_dx(k,l) * k2;
                for (int l=DIM; l--; ) H1u += this->fIntegData.fDPhiX0(i,l) * du_dx(l,k) * k2;
                Rhs[0][this->Mesh()->NLocDOF()+DIM*i+k] += (L2u + H1u) * WJ;
                
                //Lagrange multipliers residual
                double L2 = lagM_[k] * this->Mesh()->getNumericalIntegration()-> phi_(i,index) * k1;
                double H1 = 0.;
                for (int l=DIM; l--; ) H1 += this->fIntegData.fDPhiX0(i,l) * dL_dx(k,l) * k2;
                for (int l=DIM; l--; ) H1 += this->fIntegData.fDPhiX0(i,l) * dL_dx(l,k) * k2;
                Rhs[0][DIM*i+k] += (L2 + H1) * WJ;
            };
        };
    }
    if (this->Mesh()->getProblemParameters().ArlequinStab() != ENoStab){
        ArlequinStabResidual(index, this->fIntegData.fDPhiX0, this->fIntegData.fWeight, this->fIntegData.fJacA0, Rhs);
    }
}

template<class tshape>
void ElCouplingLocal<tshape>::ApplyBC(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs){
    return;
    int DIM = tshape::Dimension;

    if (this->Mesh()->getProblemParameters().ProbType() == ProblemType::EPoisson){
        for (int i = 0; i < this->Mesh()->NElNodes(); i++){
            auto connectlocal = fMeshVector[1]->ElementVec()[fLocalIndex]->getConnectivity();
            if (fMeshVector[1]->NodeVec()[connectlocal[i]] -> getConstrains(0) == 1) {
                // for (int istiff = 0; istiff < Stiffness.size(); istiff++){
                //     for (int j = 0; j < this->Mesh()->NElNodes(); j++){
                //         Stiffness[istiff](i,j) = 0.;
                //         Stiffness[istiff](j,i) = 0.;
                //     };
                // };
                Rhs[0][i] = 0.0;
                // Rhs[0][this->Mesh()->NLocDOF() + i] = 0.0;
                Rhs[1][i] = 0.0;
                // Rhs[1][this->Mesh()->NLocDOF() + i] = 0.0;
            };
        };
    } else {
        for (int i = 0; i < this->Mesh()->NElNodes(); i++){
            auto connectlocal = fMeshVector[1]->ElementVec()[fLocalIndex]->getConnectivity();
            for (int k = DIM; k--; ){
                if (fMeshVector[1]->NodeVec()[connectlocal[i]] -> getConstrains(k) == 1) {
                    // for (int istiff = 0; istiff < Stiffness.size(); istiff++){
                    //     for (int j = 0; j < this->Mesh()->NElNodes()*DIM; j++){
                    //         Stiffness[istiff](DIM*i+k,j) = 0.;
                    //         Stiffness[istiff](j,DIM*i+k) = 0.;
                    //     };
                    // };
                    Rhs[0][i] = 0.0;
                    // Rhs[0][this->Mesh()->NLocDOF() + i] = 0.0;
                    Rhs[1][i] = 0.0;
                    // Rhs[1][this->Mesh()->NLocDOF() + i] = 0.0;
                }
            };
        };
    }
}

template<class tshape>
void ElCouplingLocal<tshape>::ArlequinStabStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, std::vector<MatrixDouble> &Stiffness){

    int DIM = tshape::Dimension;
    int DEG = this->Mesh()->GetDefaultOrder();
    VecDouble xsi(DIM);

    int dimddphi = DIM == 2 ? 3 : 6;
    MatrixDouble ddphi_dx(this->Mesh()->NElNodes(),dimddphi);

    IntegQuadrature nQuad(DIM,DEG);
    auto force = this->Mesh()->getProblemParameters().getForcingFunction();

    //Defines the integration points adimentional coordinates
    for (int k = 0; k < DIM; k++) xsi[k] = nQuad.PointList(index,k);
    
    VecDouble xna_(DIM);
    VecDouble forcingF(DIM);
    for (int i = 0; i < DIM; i++) xna_[i] = this->fIntPointCoordinates(index,i);
    if (force) force(xna_,forcingF);

    //Computes the jacobian matrix
    MatrixDouble ainv_(DIM,DIM);
    // ComputeJacobian(xsi, ainv_, djac_, index);
    PanicButton();//Update it to compute in integPointData.
    this->getHighOrderSpatialDerivatives(xsi, ainv_, dphi_dx, ddphi_dx);

    double WJ = djac_ * weight_ * fMeshVector[1]->ElementVec()[fLocalIndex]->getIntegPointWeightFunction(index); 

    if (this->Mesh()->getProblemParameters().ProbType() == EPoisson){
        for (int i = 0; i < this->Mesh()->NElNodes(); i++){
            double shapeFi = this->Mesh()->getNumericalIntegration()-> phi_(i,index);
            for (int j = 0; j < this->Mesh()->NElNodes(); j++){        
                double shapeFj = this->Mesh()->getNumericalIntegration()-> phi_(j,index);
                //ARLEQUIN STABILIZATION TERMS
                double LL = 0.;
                for (int m = DIM; m--; ) LL += dphi_dx(i,m) * dphi_dx(j,m);
                // LL = shapeFi * shapeFj;

                Stiffness[1](i,j) -= 2. * (LL) * weight_ * djac_;
                
                //High order derivative
                double LH = 0.;
                // for (int m = DIM; m--; ) LH -= dphi_dx(i,m) * (ddphi_dx(j,0) + ddphi_dx(j,1));
                Stiffness[2](i,j) += LH * WJ; 
            };
        };    
    } else if (this->Mesh()->getProblemParameters().ProbType() == EElastic) {
        
        MatrixDouble matD(3,2*this->Mesh()->NElNodes());
        matD.setZero();
        MatrixDouble Hooke(3,3);
        Hooke.setZero();
        // For EPT
        double elastic_ = this->Mesh()->getProblemParameters().GetYoungModulus();
        double poisson_ = this->Mesh()->getProblemParameters().GetPoissonRatio();
        double k = elastic_ / (1. - poisson_ * poisson_);
        Hooke(0,0) = k;
        Hooke(0,1) = k * poisson_;
        Hooke(1,0) = k * poisson_;
        Hooke(1,1) = k;
        Hooke(2,2) = k * (1. - poisson_) * 0.5;
        // Hooke.setIdentity();
        for (int j = 0; j < this->Mesh()->NElNodes(); j++){
            matD(0,DIM*j  ) = dphi_dx(j,0);
            matD(1,DIM*j+1) = dphi_dx(j,1);
            matD(2,DIM*j  ) = dphi_dx(j,1);
            matD(2,DIM*j+1) = dphi_dx(j,0);
        }
        
        Stiffness[1] -= matD.transpose() * Hooke * matD * WJ;
    } else {
        PanicButton();
        // for (int i = 0; i < this->Mesh()->NElNodes(); i++){
        //     for (int j = 0; j < this->Mesh()->NElNodes(); j++){ 
        //         //ARLEQUIN STABILIZATION TERMS
        //         double LL = 0.;
        //         for (int m = DIM; m--; ) LL += dphi_dx(i,m) * dphi_dx(j,m);

        //         Stiffness[1](i,j) += (LL) * weight_ * djac_;
        //         for (int k = DIM; k--; ) Stiffness[1](DIM*i+k,DIM*j+k) += LL * weight_ * djac_;
                
        //         //High order derivative
        //         double LH = 0.;
        //         // for (int m = DIM; m--; ) LH -= dphi_dx(i,m) * (ddphi_dx(j,0) + ddphi_dx(j,1));
        //         Stiffness[2](i,j) += LH * WJ;   

        //     }
        // }     
                // PanicButton();
        //         //ARLEQUIN STABILIZATION TERMS
        //         double AM = 0.;
        //         double Lpx = 0.; double Lpy = 0.;
        //         double LC = 0.; double LL = 0.;

        //         AM = this->Mesh()->getNumericalIntegration()->phi_(i,index) * this->Mesh()->getNumericalIntegration()->phi_(j,index) * tARLQ_ * wna_* alpha_m;

        //         // LL = -2 * phi_[i] * phi_[j] * tARLQ_ / dens_;
        //         for (int m = DIM; m--; ) LL += dphi_dx(i,m) * dphi_dx(j,m) * tARLQ_ / dens_;

        //         // Lpx = -(dphi_dx[0][i] * dp_dxx + dphi_dx[1][i] * dp_dxy) * intPointWeightFunction(index)
        //         //      * tARLQ_ / dens_;
        //         // Lpy = -(dphi_dx[0][i] * dp_dxy + dphi_dx[1][i] * dp_dyy) * intPointWeightFunction(index)
        //         //      * tARLQ_ / dens_;

        //         // LC = phi_[i] * ((una_ - umesh_) * dphi_dx[0][i] + (vna_ - vmesh_) * dphi_dx[1][i]) * phi_[j] * tARLQ_ * intPointWeightFunction(index);

        //         for (int k = DIM; k--; )
        //             arlequinStab(DIM*i+k,DIM*j+k) += (AM + LL) * weight_ * djac_;
                

        //         // LC = -(dphi_dx[0][i]*(du_dx*dphi_dx[0][j] + dv_dx*dphi_dx[1][j]) +
        //         //        dphi_dx[1][i]*(du_dy*dphi_dx[0][j] + dv_dy*dphi_dx[1][j]) + 
        //         //        dphi_dx[0][i]*(u_*ddphi_dx(0,0)(j) + v_*ddphi_dx(0,1)(j)) + 
        //         //        dphi_dx[1][i]*(u_*ddphi_dx(1,0)(j) + v_*ddphi_dx(1,1)(j))) * tARLQ_ * intPointWeightFunction(index);

        //         // Lpx = 0.; Lpy = 0.;

        //         // Lpx = (dphi_dx[0][i] * ddphi_dx(0,0)(j) + 
        //         //        dphi_dx[1][i] * ddphi_dx(0,1)(j)) * tARLQ_ * intPointWeightFunction(index);
        //         // Lpy = (dphi_dx[0][i] * ddphi_dx(1,0)(j) + 
        //         //        dphi_dx[1][i] * ddphi_dx(1,1)(j)) * tARLQ_ * intPointWeightFunction(index);




        //         // laplMatrix[2*i  ][2*j  ] += (LC + AM) * weight_ * djac_;
        //         // laplMatrix[2*i+1][2*j+1] += (LC + AM) * weight_ * djac_;

        //         // laplMatrix[2*i  ][12+j] += Lpx * weight_ * djac_;
        //         // laplMatrix[2*i+1][12+j] += Lpy * weight_ * djac_;
        //         // laplMatrix[12+j][2*i  ] += Lpx * weight_ * djac_;
        //         // laplMatrix[12+j][2*i+1] += Lpy * weight_ * djac_;

        //     };

        //     //ARLEQUIN STABILIZATION TERMS
        //     double LCx = 0.; double LCy = 0.;
        //     double LPx = 0.; double LPy = 0.;

        //     for (int k = DIM; k--; ){
        //         double LLx = 0.;
        //         for (int m = DIM; m--; ) LLx -= dphi_dx(i,m) * dL_dx(k,m)/wna_ * tARLQ_ / dens_;
        //         double Amx = - this->Mesh()->getNumericalIntegration()->phi_(i,index) * am_[k] * tARLQ_;
        //         arlequinStabVector[DIM*i+k] += (Amx + LLx) * weight_ * djac_ * wna_;
        //     }

        //         // LLx = -(dphi_dx[0][i] * (dLx_dx/wna_) + dphi_dx[1][i] * (dLx_dy/wna_)) * tARLQ_ / dens_;
        //         // LLy = -(dphi_dx[0][i] * (dLy_dx/wna_) + dphi_dx[1][i] * (dLy_dy/wna_)) * tARLQ_ / dens_;
        // //}else{
        //         // LLx = -(dphi_dx[0][i] * (dLx_dx/intPointWeightFunction(index) - duna_dx*duna_dx - dvna_dx*duna_dy) + 
        //         //         dphi_dx[1][i] * (dLx_dy/intPointWeightFunction(index) - duna_dy*duna_dx - dvna_dy*duna_dy)) * tARLQ_ / dens_;
        //         // LLy = -(dphi_dx[0][i] * (dLy_dx/intPointWeightFunction(index) - duna_dx*dvna_dx - dvna_dx*dvna_dy) + 
        //         //         dphi_dx[1][i] * (dLy_dy/intPointWeightFunction(index) - duna_dy*dvna_dx - dvna_dy*dvna_dy)) * tARLQ_ / dens_;  
        //         // LLx = -(dphi_dx[0][i] * (dLx_dx/intPointWeightFunction(index) - dp_dxx) + 
        //         //         dphi_dx[1][i] * (dLx_dy/intPointWeightFunction(index) - dp_dxy)) * tARLQ_ / dens_;
        //         // LLy = -(dphi_dx[0][i] * (dLy_dx/intPointWeightFunction(index) - dp_dxy) + 
        //         //         dphi_dx[1][i] * (dLy_dy/intPointWeightFunction(index) - dp_dyy)) * tARLQ_ / dens_;  

        //         // if (iTimeStep > 10){
        //         //     LLx +=  -(dphi_dx[0][i]*dax_dx + dphi_dx[1][i]*dax_dy) * tARLQ_;
        //         //     LLx +=  -(dphi_dx[0][i]*day_dx + dphi_dx[1][i]*day_dy) * tARLQ_ ;
        //         // };

        //     //}
            
        //     // arlequinStabVector[2*i  ] += (Amx + LLx) * weight_ * djac_ * wna_;
        //     // arlequinStabVector[2*i+1] += (Amy + LLy) * weight_ * djac_ * wna_;
        // };
    }
    





};

template<class tshape>
void ElCouplingLocal<tshape>::ArlequinStabResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, std::vector<VecDouble> &Rhs){
    auto force = this->Mesh()->getProblemParameters().getForcingFunction();    
    int DIM = tshape::Dimension;
    VecDouble xna_(DIM);
    VecDouble forcingF(DIM);
    for (int i = 0; i < DIM; i++) xna_[i] = this->fIntPointCoordinates(index,i);
    if (force) force(xna_,forcingF);

    auto ellocal = fMeshVector[1]->ElementVec()[fLocalIndex];
    double WJ = djac_ * weight_ * ellocal->getIntegPointWeightFunction(index); 

    //Lagrange Multiplier
    VecDouble lagM_(this->Mesh()->NState()+1);
    this->interpolateSolution(index, lagM_);

    //Lagrange Multiplier Derivatives
    MatrixDouble dL_dx(this->Mesh()->NState()+1,DIM);
    this->interpolateSolDerivatives(dL_dx);


    if (this->Mesh()->getProblemParameters().ProbType() == EPoisson){
        for (int i = 0; i < this->Mesh()->NElNodes(); i++){
            //ARLEQUIN STABILIZATION TERMS
            double LLx = 0.;
            double LF = 0.;
            double shapeFi = this->Mesh()->getNumericalIntegration()-> phi_(i,index);
            for (int m = DIM; m--; ) LLx -= 2.*dphi_dx(i,m) * dL_dx(0,m);
            // double shapeFi = this->Mesh()->getNumericalIntegration()-> phi_(i,index);
            // for (int m = DIM; m--; ) LLx += 2*shapeFi* lagM_[0];
            // for (int m = DIM; m--; ) LF +=  dphi_dx(i,m) * forcingF[0] * fMeshVector[1]->ElementVec()[fLocalIndex]->getIntegPointWeightFunction(index);
            // double LF2 = -xna_[0]*xna_[0]*ellocal->getIntegPointWeightFunction(index)*shapeFi;
            // Rhs[1][i] += (LF2) * weight_ * djac_;
            Rhs[1][this->Mesh()->NElNodes()+i] += (LLx + LF) * weight_ * djac_;
        };    
    } else if (this->Mesh()->getProblemParameters().ProbType() == EElastic){
        VecDouble fieldForce = this->Mesh()->getProblemParameters().GetFieldForce();
        MatrixDouble matD(3,2*this->Mesh()->NElNodes());
        matD.setZero();
        MatrixDouble Hooke(3,3);
        Hooke.setZero();
        // For EPT
        double elastic_ = this->Mesh()->getProblemParameters().GetYoungModulus();
        double poisson_ = this->Mesh()->getProblemParameters().GetPoissonRatio();
        double k = elastic_ / (1. - poisson_ * poisson_);
        Hooke(0,0) = k;
        Hooke(0,1) = k * poisson_;
        Hooke(1,0) = k * poisson_;
        Hooke(1,1) = k;
        Hooke(2,2) = k * (1. - poisson_) * 0.5;
        // Hooke.setIdentity();
        for (int j = 0; j < this->Mesh()->NElNodes(); j++){
            matD(0,DIM*j  ) = dphi_dx(j,0);
            matD(1,DIM*j+1) = dphi_dx(j,1);
            matD(2,DIM*j  ) = dphi_dx(j,1);
            matD(2,DIM*j+1) = dphi_dx(j,0);
        }
        
        VecDouble strain(3);
        strain.setZero();
        strain[0] = dL_dx(0,0);
        strain[1] = dL_dx(1,1);
        strain[2] = dL_dx(0,1)+dL_dx(1,0);
        auto aux = matD.transpose() * Hooke * strain * WJ;

        VecDouble forcingF(2);
        VecDouble x_ = this->getIntegPointCoordinatesValue(index);
        if (force) force(x_,forcingF);

        for (int i = this->Mesh()->NElNodes(); i--; ){
            Rhs[1][2*this->Mesh()->NElNodes()+2*i  ] += aux[2*i];
            Rhs[1][2*this->Mesh()->NElNodes()+2*i+1] += aux[2*i+1];
            double shapeFi = this->Mesh()->getNumericalIntegration()-> phi_(i,index);
            //External force
            double Fx = (fieldForce[0] + forcingF[0]) * shapeFi;
            double Fy = (fieldForce[1] + forcingF[1]) * shapeFi;
            Rhs[1][2*this->Mesh()->NElNodes()+2*i  ] += Fx * WJ;
            Rhs[1][2*this->Mesh()->NElNodes()+2*i+1] += Fy * WJ;
        };

    } else {
        return;
        PanicButton();
    }
};



#include "ShapeHexahedron.h"
#include "ShapeOneD.h"
#include "ShapeQuadrilateral.h"
#include "ShapePoint.h"
#include "ShapeTetrahedron.h"
#include "ShapeTriangle.h"

template class ElCouplingLocal<ShapePoint>;
template class ElCouplingLocal<ShapeOneD>;
template class ElCouplingLocal<ShapeTriangle>;
template class ElCouplingLocal<ShapeQuadrilateral>;
template class ElCouplingLocal<ShapeTetrahedron>;
template class ElCouplingLocal<ShapeHexahedron>;