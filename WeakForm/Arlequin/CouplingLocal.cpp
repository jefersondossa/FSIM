#include "CouplingLocal.h"

CouplingLocal::CouplingLocal(int dim, int64_t fineindex, CompMesh* meshlocal, double k0, double k1) : WeakForm(){
    fLocalMesh = meshlocal;
    fLocalIndex = fineindex;
    fDimension = dim;
    fK0 = k0;
    fK1 = k1;
    this->fType = WeakFormType::kCouplingLocal;
};

 
void CouplingLocal::ComputeStiffness(int &index, IntPointData &data, std::vector<MatrixDouble> &Stiffness){

    //Stiffness[0] is always the L2-H1 contribution;
    //Stiffness[1] is always the diagonal contribution to the global stiffness matrix
    //The remaining positions corresponds to the other contributions, 
    //which will be added to stiffness[0] in the global matrix and depends on the problem.
    if (!data.fNeedsDSol || !data.fNeedsSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
        data.fNeedsSol = true;
        data.fSol.resize(fNState);
    }

    int DIM = fDimension;
    int nphi = data.fPhi.size();
    double WJ = data.fWeight * data.fJacA0;
    int nstate = NState();
    // if (this->Mesh()->getProblemParameters().ProbType() == ProblemType::EPoisson){
    
    for (int i = 0; i < nphi; i++){
        for (int j = 0; j < nphi; j++){
            double l2 = data.fPhi[i] * data.fPhi[j] * WJ * fK0;
                for (int istate = 0; istate < nstate; istate++){
                // L2 COUPLING OPERATOR
                Stiffness[0](nstate*i+istate,nstate*j+istate) -= l2;

                // H1 COUPLING OPERATOR
                if (fK1 > 0){
                    for (int k = DIM; k--;  ){
                        for (int l = DIM; l--; ){
                            double K = data.fDPhi(i,l) * data.fDPhi(j,k) * fK1;
                            if (k==l) for (int m = DIM; m--; ) K += data.fDPhiX0(m,i) * data.fDPhiX0(m,j);
                            Stiffness[0](DIM*i+k,DIM*j+l) -= K * WJ;
                        }
                    }
                }                
 
                // for (int l = 0; l < DIM; l++){
                //     
                //     double K = data.fDPhiX0(i,l) * data.fDPhiX0(j,l);
                //     Stiffness[0](i,j) -= K * WJ * fK1;
                // };
            };
        };
    }
    // } else {
    //     for (int i = 0; i < nphi; i++){
    //         for (int j = 0; j < nphi; j++){
    //             double l2 = data.fPhi[i] * data.fPhi[j] * WJ * fK0;
    //             for (int k = 0; k < DIM; k++){
    //                 // L2 COUPLING OPERATOR
    //                 Stiffness[0](DIM*i+k,DIM*j+k) -= l2;
    //                 for (int l = 0; l < DIM; l++){
    //                     //H1 COUPLING OPERATOR
    //                     double K = data.fDPhiX0(i,l) * data.fDPhiX0(j,k);
    //                     if (k==l) for (int m = DIM; m--; ) K += data.fDPhiX0(i,m) * data.fDPhiX0(j,m);

    //                     Stiffness[0](DIM*i+k,DIM*j+l) -= K * WJ * fK1;
    //                 };//l
    //             };//k
    //         };//j
    //     };//i
    // };//if Poisson

    // if (this->Mesh()->getProblemParameters().ArlequinStab() != ENoStab){
    //     ArlequinStabStiffness(index, data.fDPhiX0, data.fWeight, data.fJacA0, Stiffness);
    // }

};

 
void CouplingLocal::ComputeResidual(int &index, IntPointData &data, std::vector<VecDouble> &Rhs){


    int DIM = fDimension;
    int nphi = data.fPhi.size();
    auto *ellocal = fLocalMesh->ElementVec()[fLocalIndex];

    //Velocity
    VecDouble u_(fNState+1);
    ellocal->interpolateSolution(data.fPhi, u_);

    //Velocity Derivatives
    MatrixDouble du_dx(fNState+1,DIM);
    ellocal->interpolateSolDerivatives(data.fDPhiX0,du_dx);
    
    double WJ = data.fWeight * data.fJacA0;

    // if (this->Mesh()->getProblemParameters().ProbType() == ProblemType::EPoisson){
        for (int i = 0; i < nphi; i++){
            //Solution residual
            double L2u = u_[0] * data.fPhi[i] * fK0;
            double H1u = 0.;
            for (int l=DIM; l--; ) H1u += data.fDPhiX0(l,i) * du_dx(0,l) * fK1;
            Rhs[0][nphi*fNState+i] += (L2u + H1u) * WJ;
            
            // Lagrange multipliers residual
            double L2 = data.fSol[0] * data.fPhi[i] * fK0;
            double H1 = 0.;
            for (int l=DIM; l--; ) H1 += data.fDPhiX0(l,i) * data.fDSolDx(0,l) * fK1;
            Rhs[0][i] += (L2 + H1) * WJ;
        };
    // } else {
    //     for (int i = 0; i < nphi; i++){
    //         for (int k = 0; k < DIM; k++){
    //             //Solution residual
    //             double L2u = u_[k] * data.fPhi[i] * fK0;
    //             double H1u = 0.;
    //             for (int l=DIM; l--; ) H1u += data.fDPhiX0(i,l) * du_dx(k,l) * fK1;
    //             for (int l=DIM; l--; ) H1u += data.fDPhiX0(i,l) * du_dx(l,k) * fK1;
    //             Rhs[0][Element::NLocDOF()+DIM*i+k] += (L2u + H1u) * WJ;
                
    //             //Lagrange multipliers residual
    //             double L2 = lagM_[k] * data.fPhi[i] * fK0;
    //             double H1 = 0.;
    //             for (int l=DIM; l--; ) H1 += data.fDPhiX0(i,l) * dL_dx(k,l) * fK1;
    //             for (int l=DIM; l--; ) H1 += data.fDPhiX0(i,l) * dL_dx(l,k) * fK1;
    //             Rhs[0][DIM*i+k] += (L2 + H1) * WJ;
    //         };
    //     };
    // }
    // if (this->Mesh()->getProblemParameters().ArlequinStab() != ENoStab){
    //     ArlequinStabResidual(index, data.fDPhiX0, data.fWeight, data.fJacA0, Rhs);
    // }
}

 
// void CouplingLocal::ArlequinStabStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, std::vector<MatrixDouble> &Stiffness){

//     int DIM = fDimension;
//     VecDouble xsi(DIM);

//     int dimddphi = DIM == 2 ? 3 : 6;
//     MatrixDouble ddphi_dx(nphi,dimddphi);

//     auto force = fForceFunction;

//     //Defines the integration points adimentional coordinates
//     for (int k = 0; k < DIM; k++) xsi[k] = data.fAdimCoord;
    
//     VecDouble xna_(DIM);
//     VecDouble forcingF(DIM);
//     for (int i = 0; i < DIM; i++) xna_[i] = this->fIntPointCoordinates(index,i);
//     if (force) force(xna_,forcingF);

//     //Computes the jacobian matrix
//     MatrixDouble ainv_(DIM,DIM);
//     // ComputeJacobian(xsi, ainv_, djac_, index);
//     PanicButton();//Update it to compute in integPointData.
//     // this->getHighOrderSpatialDerivatives(xsi, ainv_, dphi_dx, ddphi_dx);

//     double WJ = djac_ * weight_ * fMeshVector[1]->ElementVec()[fLocalIndex]->getIntegPointWeightFunction(index); 

//     if (this->Mesh()->getProblemParameters().ProbType() == EPoisson){
//         for (int i = 0; i < nphi; i++){
//             double shapeFi = data.fPhi[i];
//             for (int j = 0; j < nphi; j++){        
//                 double shapeFj = data.fPhi[j];
//                 //ARLEQUIN STABILIZATION TERMS
//                 double LL = 0.;
//                 for (int m = DIM; m--; ) LL += dphi_dx(i,m) * dphi_dx(j,m);
//                 // LL = shapeFi * shapeFj;

//                 Stiffness[1](i,j) -= 2. * (LL) * weight_ * djac_;
                
//                 //High order derivative
//                 double LH = 0.;
//                 // for (int m = DIM; m--; ) LH -= dphi_dx(i,m) * (ddphi_dx(j,0) + ddphi_dx(j,1));
//                 Stiffness[2](i,j) += LH * WJ; 
//             };
//         };    
//     } else if (this->Mesh()->getProblemParameters().ProbType() == EElastic) {
        
//         MatrixDouble matD(3,2*nphi);
//         matD.setZero();
//         MatrixDouble Hooke(3,3);
//         Hooke.setZero();
//         // For EPT
//         double elastic_ = this->Mesh()->getProblemParameters().GetYoungModulus();
//         double poisson_ = this->Mesh()->getProblemParameters().GetPoissonRatio();
//         double k = elastic_ / (1. - poisson_ * poisson_);
//         Hooke(0,0) = k;
//         Hooke(0,1) = k * poisson_;
//         Hooke(1,0) = k * poisson_;
//         Hooke(1,1) = k;
//         Hooke(2,2) = k * (1. - poisson_) * 0.5;
//         // Hooke.setIdentity();
//         for (int j = 0; j < nphi; j++){
//             matD(0,DIM*j  ) = dphi_dx(j,0);
//             matD(1,DIM*j+1) = dphi_dx(j,1);
//             matD(2,DIM*j  ) = dphi_dx(j,1);
//             matD(2,DIM*j+1) = dphi_dx(j,0);
//         }
        
//         Stiffness[1] -= matD.transpose() * Hooke * matD * WJ;
//     } else {
//         PanicButton();
//         // for (int i = 0; i < tshape::NElNodes; i++){
//         //     for (int j = 0; j < tshape::NElNodes; j++){ 
//         //         //ARLEQUIN STABILIZATION TERMS
//         //         double LL = 0.;
//         //         for (int m = DIM; m--; ) LL += dphi_dx(i,m) * dphi_dx(j,m);

//         //         Stiffness[1](i,j) += (LL) * weight_ * djac_;
//         //         for (int k = DIM; k--; ) Stiffness[1](DIM*i+k,DIM*j+k) += LL * weight_ * djac_;
                
//         //         //High order derivative
//         //         double LH = 0.;
//         //         // for (int m = DIM; m--; ) LH -= dphi_dx(i,m) * (ddphi_dx(j,0) + ddphi_dx(j,1));
//         //         Stiffness[2](i,j) += LH * WJ;   

//         //     }
//         // }     
//                 // PanicButton();
//         //         //ARLEQUIN STABILIZATION TERMS
//         //         double AM = 0.;
//         //         double Lpx = 0.; double Lpy = 0.;
//         //         double LC = 0.; double LL = 0.;

//         //         AM = data.fPhi[i] * data.fPhi[j] * tARLQ_ * wna_* alpha_m;

//         //         // LL = -2 * phi_[i] * phi_[j] * tARLQ_ / dens_;
//         //         for (int m = DIM; m--; ) LL += dphi_dx(i,m) * dphi_dx(j,m) * tARLQ_ / dens_;

//         //         // Lpx = -(dphi_dx[0][i] * dp_dxx + dphi_dx[1][i] * dp_dxy) * data.fWeightFunction(index)
//         //         //      * tARLQ_ / dens_;
//         //         // Lpy = -(dphi_dx[0][i] * dp_dxy + dphi_dx[1][i] * dp_dyy) * data.fWeightFunction(index)
//         //         //      * tARLQ_ / dens_;

//         //         // LC = phi_[i] * ((una_ - umesh_) * dphi_dx[0][i] + (vna_ - vmesh_) * dphi_dx[1][i]) * phi_[j] * tARLQ_ * data.fWeightFunction(index);

//         //         for (int k = DIM; k--; )
//         //             arlequinStab(DIM*i+k,DIM*j+k) += (AM + LL) * weight_ * djac_;
                

//         //         // LC = -(dphi_dx[0][i]*(du_dx*dphi_dx[0][j] + dv_dx*dphi_dx[1][j]) +
//         //         //        dphi_dx[1][i]*(du_dy*dphi_dx[0][j] + dv_dy*dphi_dx[1][j]) + 
//         //         //        dphi_dx[0][i]*(u_*ddphi_dx(0,0)(j) + v_*ddphi_dx(0,1)(j)) + 
//         //         //        dphi_dx[1][i]*(u_*ddphi_dx(1,0)(j) + v_*ddphi_dx(1,1)(j))) * tARLQ_ * data.fWeightFunction(index);

//         //         // Lpx = 0.; Lpy = 0.;

//         //         // Lpx = (dphi_dx[0][i] * ddphi_dx(0,0)(j) + 
//         //         //        dphi_dx[1][i] * ddphi_dx(0,1)(j)) * tARLQ_ * data.fWeightFunction(index);
//         //         // Lpy = (dphi_dx[0][i] * ddphi_dx(1,0)(j) + 
//         //         //        dphi_dx[1][i] * ddphi_dx(1,1)(j)) * tARLQ_ * data.fWeightFunction(index);




//         //         // laplMatrix[2*i  ][2*j  ] += (LC + AM) * weight_ * djac_;
//         //         // laplMatrix[2*i+1][2*j+1] += (LC + AM) * weight_ * djac_;

//         //         // laplMatrix[2*i  ][12+j] += Lpx * weight_ * djac_;
//         //         // laplMatrix[2*i+1][12+j] += Lpy * weight_ * djac_;
//         //         // laplMatrix[12+j][2*i  ] += Lpx * weight_ * djac_;
//         //         // laplMatrix[12+j][2*i+1] += Lpy * weight_ * djac_;

//         //     };

//         //     //ARLEQUIN STABILIZATION TERMS
//         //     double LCx = 0.; double LCy = 0.;
//         //     double LPx = 0.; double LPy = 0.;

//         //     for (int k = DIM; k--; ){
//         //         double LLx = 0.;
//         //         for (int m = DIM; m--; ) LLx -= dphi_dx(i,m) * dL_dx(k,m)/wna_ * tARLQ_ / dens_;
//         //         double Amx = - data.fPhi[i] * am_[k] * tARLQ_;
//         //         arlequinStabVector[DIM*i+k] += (Amx + LLx) * weight_ * djac_ * wna_;
//         //     }

//         //         // LLx = -(dphi_dx[0][i] * (dLx_dx/wna_) + dphi_dx[1][i] * (dLx_dy/wna_)) * tARLQ_ / dens_;
//         //         // LLy = -(dphi_dx[0][i] * (dLy_dx/wna_) + dphi_dx[1][i] * (dLy_dy/wna_)) * tARLQ_ / dens_;
//         // //}else{
//         //         // LLx = -(dphi_dx[0][i] * (dLx_dx/data.fWeightFunction(index) - duna_dx*duna_dx - dvna_dx*duna_dy) + 
//         //         //         dphi_dx[1][i] * (dLx_dy/data.fWeightFunction(index) - duna_dy*duna_dx - dvna_dy*duna_dy)) * tARLQ_ / dens_;
//         //         // LLy = -(dphi_dx[0][i] * (dLy_dx/data.fWeightFunction(index) - duna_dx*dvna_dx - dvna_dx*dvna_dy) + 
//         //         //         dphi_dx[1][i] * (dLy_dy/data.fWeightFunction(index) - duna_dy*dvna_dx - dvna_dy*dvna_dy)) * tARLQ_ / dens_;  
//         //         // LLx = -(dphi_dx[0][i] * (dLx_dx/data.fWeightFunction(index) - dp_dxx) + 
//         //         //         dphi_dx[1][i] * (dLx_dy/data.fWeightFunction(index) - dp_dxy)) * tARLQ_ / dens_;
//         //         // LLy = -(dphi_dx[0][i] * (dLy_dx/data.fWeightFunction(index) - dp_dxy) + 
//         //         //         dphi_dx[1][i] * (dLy_dy/data.fWeightFunction(index) - dp_dyy)) * tARLQ_ / dens_;  

//         //         // if (iTimeStep > 10){
//         //         //     LLx +=  -(dphi_dx[0][i]*dax_dx + dphi_dx[1][i]*dax_dy) * tARLQ_;
//         //         //     LLx +=  -(dphi_dx[0][i]*day_dx + dphi_dx[1][i]*day_dy) * tARLQ_ ;
//         //         // };

//         //     //}
            
//         //     // arlequinStabVector[2*i  ] += (Amx + LLx) * weight_ * djac_ * wna_;
//         //     // arlequinStabVector[2*i+1] += (Amy + LLy) * weight_ * djac_ * wna_;
//         // };
//     }
    





// };

 
// void CouplingLocal::ArlequinStabResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, std::vector<VecDouble> &Rhs){
//     auto force = this->Mesh()->getProblemParameters().getForcingFunction();    
//     int DIM = tshape::Dimension;
//     VecDouble xna_(DIM);
//     VecDouble forcingF(DIM);
//     for (int i = 0; i < DIM; i++) xna_[i] = this->fIntPointCoordinates(index,i);
//     if (force) force(xna_,forcingF);

//     auto ellocal = fMeshVector[1]->ElementVec()[fLocalIndex];
//     double WJ = djac_ * weight_ * ellocal->getIntegPointWeightFunction(index); 

//     //Lagrange Multiplier
//     VecDouble lagM_(this->Mesh()->NState()+1);
//     this->interpolateSolution(index, lagM_);

//     //Lagrange Multiplier Derivatives
//     MatrixDouble dL_dx(this->Mesh()->NState()+1,DIM);
//     this->interpolateSolDerivatives(dL_dx);


//     if (this->Mesh()->getProblemParameters().ProbType() == EPoisson){
//         for (int i = 0; i < tshape::NElNodes; i++){
//             //ARLEQUIN STABILIZATION TERMS
//             double LLx = 0.;
//             double LF = 0.;
//             double shapeFi = data.fPhi[i];
//             for (int m = DIM; m--; ) LLx -= 2.*dphi_dx(i,m) * dL_dx(0,m);
//             // double shapeFi = data.fPhi[i];
//             // for (int m = DIM; m--; ) LLx += 2*shapeFi* lagM_[0];
//             // for (int m = DIM; m--; ) LF +=  dphi_dx(i,m) * forcingF[0] * fMeshVector[1]->ElementVec()[fLocalIndex]->getIntegPointWeightFunction(index);
//             // double LF2 = -xna_[0]*xna_[0]*ellocal->getIntegPointWeightFunction(index)*shapeFi;
//             // Rhs[1][i] += (LF2) * weight_ * djac_;
//             Rhs[1][tshape::NElNodes+i] += (LLx + LF) * weight_ * djac_;
//         };    
//     } else if (this->Mesh()->getProblemParameters().ProbType() == EElastic){
//         VecDouble fieldForce = this->Mesh()->getProblemParameters().GetFieldForce();
//         MatrixDouble matD(3,2*tshape::NElNodes);
//         matD.setZero();
//         MatrixDouble Hooke(3,3);
//         Hooke.setZero();
//         // For EPT
//         double elastic_ = this->Mesh()->getProblemParameters().GetYoungModulus();
//         double poisson_ = this->Mesh()->getProblemParameters().GetPoissonRatio();
//         double k = elastic_ / (1. - poisson_ * poisson_);
//         Hooke(0,0) = k;
//         Hooke(0,1) = k * poisson_;
//         Hooke(1,0) = k * poisson_;
//         Hooke(1,1) = k;
//         Hooke(2,2) = k * (1. - poisson_) * 0.5;
//         // Hooke.setIdentity();
//         for (int j = 0; j < tshape::NElNodes; j++){
//             matD(0,DIM*j  ) = dphi_dx(j,0);
//             matD(1,DIM*j+1) = dphi_dx(j,1);
//             matD(2,DIM*j  ) = dphi_dx(j,1);
//             matD(2,DIM*j+1) = dphi_dx(j,0);
//         }
        
//         VecDouble strain(3);
//         strain.setZero();
//         strain[0] = dL_dx(0,0);
//         strain[1] = dL_dx(1,1);
//         strain[2] = dL_dx(0,1)+dL_dx(1,0);
//         auto aux = matD.transpose() * Hooke * strain * WJ;

//         VecDouble forcingF(2);
//         VecDouble x_ = this->getIntegPointCoordinatesValue(index);
//         if (force) force(x_,forcingF);

//         for (int i = nphi; i--; ){
//             Rhs[1][2*nphi+2*i  ] += aux[2*i];
//             Rhs[1][2*nphi+2*i+1] += aux[2*i+1];
//             double shapeFi = data.fPhi[i];
//             //External force
//             double Fx = (fieldForce[0] + forcingF[0]) * shapeFi;
//             double Fy = (fieldForce[1] + forcingF[1]) * shapeFi;
//             Rhs[1][2*nphi+2*i  ] += Fx * WJ;
//             Rhs[1][2*nphi+2*i+1] += Fy * WJ;
//         };

//     } else {
//         return;
//         PanicButton();
//     }
// };
