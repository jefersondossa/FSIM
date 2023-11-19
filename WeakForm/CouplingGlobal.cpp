#include "CouplingGlobal.h"


void CouplingGlobal::ComputeStiffness(int &index, IntPointData &data, std::vector<MatrixDouble> &Stiffness){
    
    //Stiffness[0] is always the L2-H1 contribution;
    //Stiffness[1] is always the diagonal contribution to the global stiffness matrix
    //The remaining positions corresponds to the other contributions, 
    //which will be added to stiffness[0] in the global matrix and depends on the problem.

    if ((fGlobalElemCorresp[index] != fGlobalIndex)) return;// This contribution belongs to another global element
    
    if (!data.fNeedsDSol || !data.fNeedsSol){
        data.fNeedsDSol = true;
        data.fDSolDx.resize(fNState,fDimension);
        data.fNeedsSol = true;
        data.fSol.resize(fNState);
    }

    int nphi = data.fPhi.size();
    int DIM = fDimension;

    double WJ = data.fWeight * data.fJacA0;
    
    VecDouble XsiGlobal(DIM);
    for (int k = 0; k < DIM; k++) XsiGlobal[k] = fGlobalXsi(index,k);
    //Computes the coarse mesh shape functions
    auto *elglobal = fGlobalMesh->ElementVec()[fGlobalIndex];
    auto &dataglobal = elglobal->IntegrationData();
    dataglobal.fAdimCoord = XsiGlobal;
    // VecDouble phiGlobal(nphi);
    // tshape::Shape(XsiGlobal,phiGlobal);
    //Computes coarse mesh derivatives
    elglobal->ComputeJacobian();
    elglobal->ComputeSpatialDerivatives();

    auto dphi_dx = data.fDPhiX0;
    auto dphi_dxGlobal = dataglobal.fDPhiX0;
    
    // if (this->Mesh()->getProblemParameters().ProbType() == ProblemType::EPoisson){
        for (int i = 0; i < nphi; i++){
            for (int j = 0; j < nphi; j++){
                double l2 = data.fPhi[i] * dataglobal.fPhi[j] * WJ * fK0;
                // L2 COUPLING OPERATOR
                Stiffness[0](i,j) += l2;
                for (int l = 0; l < DIM; l++){
                    //H1 COUPLING OPERATOR
                    double K = dphi_dx(i,l) * dphi_dxGlobal(j,l);
                    Stiffness[0](i,j) += K * WJ * fK1;
                };
            };
        };
    // } else {
    //     for (int i = 0; i < nphi; i++){
    //         for (int j = 0; j < nphi; j++){
    //             double l2 = data.fPhi[i] * phiGlobal[j] * WJ * fK0;
    //             for (int k = 0; k < DIM; k++){
    //                 // L2 COUPLING OPERATOR
    //                 Stiffness[0](DIM*i+k,DIM*j+k) += l2;
    //                 for (int l = 0; l < DIM; l++){
    //                     //H1 COUPLING OPERATOR
    //                     double K = dphi_dx(i,l) * dphi_dx(j,k);
    //                     if (k==l) for (int m = DIM; m--; ) K += dphi_dx(i,m) * dphi_dxGlobal(j,m);

    //                     Stiffness[0](DIM*i+k,DIM*j+l) += K * WJ * fK1;
    //                 };//l
    //             };//k
    //         };//j
    //     };//i
    // };//if Poisson

    // if (this->Mesh()->getProblemParameters().ArlequinStab() != ENoStab){
    //     ArlequinStabStiffness(index, dphi_dx, phiGlobal, dphi_dxGlobal, data.fWeight, data.fJacA0, Stiffness);
    // }

};



void CouplingGlobal::ComputeResidual(int &index, IntPointData &data, std::vector<VecDouble> &Rhs){

    int DIM = fDimension;
    int nphi = data.fPhi.size();
    auto *elglobal = fGlobalMesh->ElementVec()[fGlobalIndex];

    auto &dataglobal = elglobal->IntegrationData();
    //In this case, we need the global solution
    //Velocity
    VecDouble u_(fNState+1);
    elglobal->interpolateSolution(dataglobal.fPhi, u_);

    //Velocity Derivatives
    MatrixDouble du_dx(fNState+1,DIM);
    elglobal->interpolateSolDerivatives(du_dx);

    
    auto dphi_dxGlobal = dataglobal.fDPhiX0;

    double WJ = data.fWeight * data.fJacA0;
    auto dphi_dx = data.fDPhiX0;

    // if (this->Mesh()->getProblemParameters().ProbType() == ProblemType::EPoisson){
        for (int i = 0; i < nphi; i++){
            //Solution residual
            double L2u = u_[0] * data.fPhi[i] * fK0;
            double H1u = 0.;
            for (int l=DIM; l--; ) H1u += dphi_dx(i,l) * du_dx(0,l) * fK1;
            Rhs[0][nphi*fNState+i] -= (L2u + H1u) * WJ;
            
            // Lagrange multipliers residual
            double L2 = data.fSol[0] * dataglobal.fPhi[i] * fK0;
            double H1 = 0.;
            for (int l=DIM; l--; ) H1 += dphi_dxGlobal(i,l) * data.fDSolDx(0,l) * fK1;
            Rhs[0][i] -= (L2 + H1) * WJ;
        };
    // } else {
    //     for (int i = 0; i < nphi; i++){
    //         for (int k = 0; k < DIM; k++){
    //             //Solution residual
    //             double L2u = u_[k] * data.fPhi[i] * fK0;
    //             double H1u = 0.;
    //             for (int l=DIM; l--; ) H1u += dphi_dx(i,l) * du_dx(k,l) * fK1;
    //             for (int l=DIM; l--; ) H1u += dphi_dx(i,l) * du_dx(l,k) * fK1;
    //             Rhs[0][Element::NLocDOF()+DIM*i+k] -= (L2u + H1u) * WJ;
                
    //             //Lagrange multipliers residual
    //             double L2 = lagM_[k] * phiGlobal[i] * fK0;
    //             double H1 = 0.;
    //             for (int l=DIM; l--; ) H1 += dphi_dxGlobal(i,l) * dL_dx(k,l) * fK1;
    //             for (int l=DIM; l--; ) H1 += dphi_dxGlobal(i,l) * dL_dx(l,k) * fK1;
    //             Rhs[0][DIM*i+k] -= (L2 + H1) * WJ;
    //         };
    //     };
    // }
    // if (this->Mesh()->getProblemParameters().ArlequinStab() != ENoStab){
    //     ArlequinStabResidual(index, dphi_dx, phiGlobal, dphi_dxGlobal, data.fWeight, data.fJacA0, Rhs);
    // }
}


void CouplingGlobal::ArlequinStabStiffness(int &index, MatrixDouble &dphi_dx, VecDouble &phiGlobal, MatrixDouble &dphi_dxGlobal, double &weight_, double &djac_, std::vector<MatrixDouble> &Stiffness){

    // int DIM = tshape::Dimension;
    // int DEG = this->Mesh()->GetDefaultOrder();
    // VecDouble xsi(DIM);
    // int dimddphi = DIM == 2 ? 3 : 6;
    // MatrixDouble ddphi_dxGlobal(nphi,dimddphi);

    // IntegQuadrature nQuad(DIM,DEG);
    // auto force = this->Mesh()->getProblemParameters().getForcingFunction();

    // //Defines the integration points adimentional coordinates
    // for (int k = 0; k < DIM; k++) xsi[k] = nQuad.PointList(index,k);

    // //Computes the jacobian matrix
    // double djacG_;
    // MatrixDouble ainvG_(DIM,DIM);
    // // ComputeJacobian(xsi, ainvG_, djacG_, index);
    // // getHighOrderSpatialDerivatives(xsi, ainvG_, dphi_dxGlobal, ddphi_dxGlobal);

    // double WJ = djac_ * weight_ * fMeshVector[0]->ElementVec()[fGlobalIndex]->getIntegPointWeightFunction(index); 

    // if (this->Mesh()->getProblemParameters().ProbType() == EPoisson){
    //     for (int i = 0; i < nphi; i++){
    //         for (int j = 0; j < nphi; j++){        

    //             //ARLEQUIN STABILIZATION TERMS
    //             double LL = 0.;
    //             // for (int m = DIM; m--; ) LL += dphi_dx(i,m) * dphi_dxGlobal(j,m);

    //             // Stiffness[1](i,j) += (LL) * weight_ * djac_;
                
    //             //High order derivative
    //             double LH = 0.;
    //             // for (int m = DIM; m--; ) LH -= dphi_dx(i,m) * (ddphi_dxGlobal(j,0) + ddphi_dxGlobal(j,1));
    //             Stiffness[2](i,j) += LH * WJ; 
    //         };
    //     };    
    // } else if (this->Mesh()->getProblemParameters().ProbType() == EElastic) {
    //     MatrixDouble matD(3,2*nphi);
    //     MatrixDouble matDGlob(3,2*nphi);
    //     matD.setZero();
    //     matDGlob.setZero();
    //     MatrixDouble Hooke(3,3);
    //     Hooke.setZero();
    //     // For EPT
    //     double elastic_ = this->Mesh()->getProblemParameters().GetYoungModulus();
    //     double poisson_ = this->Mesh()->getProblemParameters().GetPoissonRatio();
    //     double k = elastic_ / (1. - poisson_ * poisson_);
    //     Hooke(0,0) = k;
    //     Hooke(0,1) = k * poisson_;
    //     Hooke(1,0) = k * poisson_;
    //     Hooke(1,1) = k;
    //     Hooke(2,2) = k * (1. - poisson_) * 0.5;
    //     // Hooke.setIdentity();
    //     for (int j = 0; j < nphi; j++){
    //         matD(0,DIM*j  ) = dphi_dx(j,0);
    //         matD(1,DIM*j+1) = dphi_dx(j,1);
    //         matD(2,DIM*j  ) = dphi_dx(j,1);
    //         matD(2,DIM*j+1) = dphi_dx(j,0);
    //         matDGlob(0,DIM*j  ) = dphi_dxGlobal(j,0);
    //         matDGlob(1,DIM*j+1) = dphi_dxGlobal(j,1);
    //         matDGlob(2,DIM*j  ) = dphi_dxGlobal(j,1);
    //         matDGlob(2,DIM*j+1) = dphi_dxGlobal(j,0);
    //     }
        
    //     Stiffness[1] -= matD.transpose() * Hooke * matDGlob * WJ;

    // } else {
    //     PanicButton();
    //     // for (int i = 0; i < nphi; i++){
    //     //     for (int j = 0; j < nphi; j++){ 
    //     //         //ARLEQUIN STABILIZATION TERMS
    //     //         double LL = 0.;
    //     //         for (int m = DIM; m--; ) LL += dphi_dx(i,m) * dphi_dxGlobal(j,m);

    //     //         Stiffness[1](i,j) += (LL) * weight_ * djac_;
    //     //         for (int k = DIM; k--; ) Stiffness[1](DIM*i+k,DIM*j+k) += LL * weight_ * djac_;
                
    //     //         //High order derivative
    //     //         double LH = 0.;
    //     //         // for (int m = DIM; m--; ) LH -= dphi_dx(i,m) * (ddphi_dx(j,0) + ddphi_dx(j,1));
    //     //         Stiffness[2](i,j) += LH * WJ;   

    //     //     }
    //     // }     
    //     // for (int i = 0; i < nphi; i++){
    //         // for (int j = 0; j < nphi; j++){        
    //             // PanicButton();
    //     //         //ARLEQUIN STABILIZATION TERMS
    //     //         double AM = 0.;
    //     //         double Lpx = 0.; double Lpy = 0.;
    //     //         double LC = 0.; double LL = 0.;

    //     //         AM = data.fPhi[i] * data.fPhi[j] * tARLQ_ * wna_* alpha_m;

    //     //         // LL = -2 * phi_[i] * phi_[j] * tARLQ_ / dens_;
    //     //         for (int m = DIM; m--; ) LL += dphi_dx(i,m) * dphi_dx(j,m) * tARLQ_ / dens_;

    //     //         // Lpx = -(dphi_dx[0][i] * dp_dxx + dphi_dx[1][i] * dp_dxy) * data.fWeightFunction(index)
    //     //         //      * tARLQ_ / dens_;
    //     //         // Lpy = -(dphi_dx[0][i] * dp_dxy + dphi_dx[1][i] * dp_dyy) * data.fWeightFunction(index)
    //     //         //      * tARLQ_ / dens_;

    //     //         // LC = phi_[i] * ((una_ - umesh_) * dphi_dx[0][i] + (vna_ - vmesh_) * dphi_dx[1][i]) * phi_[j] * tARLQ_ * data.fWeightFunction(index);

    //     //         for (int k = DIM; k--; )
    //     //             arlequinStab(DIM*i+k,DIM*j+k) += (AM + LL) * weight_ * djac_;
                

    //     //         // LC = -(dphi_dx[0][i]*(du_dx*dphi_dx[0][j] + dv_dx*dphi_dx[1][j]) +
    //     //         //        dphi_dx[1][i]*(du_dy*dphi_dx[0][j] + dv_dy*dphi_dx[1][j]) + 
    //     //         //        dphi_dx[0][i]*(u_*ddphi_dx(0,0)(j) + v_*ddphi_dx(0,1)(j)) + 
    //     //         //        dphi_dx[1][i]*(u_*ddphi_dx(1,0)(j) + v_*ddphi_dx(1,1)(j))) * tARLQ_ * data.fWeightFunction(index);

    //     //         // Lpx = 0.; Lpy = 0.;

    //     //         // Lpx = (dphi_dx[0][i] * ddphi_dx(0,0)(j) + 
    //     //         //        dphi_dx[1][i] * ddphi_dx(0,1)(j)) * tARLQ_ * data.fWeightFunction(index);
    //     //         // Lpy = (dphi_dx[0][i] * ddphi_dx(1,0)(j) + 
    //     //         //        dphi_dx[1][i] * ddphi_dx(1,1)(j)) * tARLQ_ * data.fWeightFunction(index);




    //     //         // laplMatrix[2*i  ][2*j  ] += (LC + AM) * weight_ * djac_;
    //     //         // laplMatrix[2*i+1][2*j+1] += (LC + AM) * weight_ * djac_;

    //     //         // laplMatrix[2*i  ][12+j] += Lpx * weight_ * djac_;
    //     //         // laplMatrix[2*i+1][12+j] += Lpy * weight_ * djac_;
    //     //         // laplMatrix[12+j][2*i  ] += Lpx * weight_ * djac_;
    //     //         // laplMatrix[12+j][2*i+1] += Lpy * weight_ * djac_;

    //     //     };

    //     //     //ARLEQUIN STABILIZATION TERMS
    //     //     double LCx = 0.; double LCy = 0.;
    //     //     double LPx = 0.; double LPy = 0.;

    //     //     for (int k = DIM; k--; ){
    //     //         double LLx = 0.;
    //     //         for (int m = DIM; m--; ) LLx -= dphi_dx(i,m) * dL_dx(k,m)/wna_ * tARLQ_ / dens_;
    //     //         double Amx = - data.fPhi[i] * am_[k] * tARLQ_;
    //     //         arlequinStabVector[DIM*i+k] += (Amx + LLx) * weight_ * djac_ * wna_;
    //     //     }

    //     //         // LLx = -(dphi_dx[0][i] * (dLx_dx/wna_) + dphi_dx[1][i] * (dLx_dy/wna_)) * tARLQ_ / dens_;
    //     //         // LLy = -(dphi_dx[0][i] * (dLy_dx/wna_) + dphi_dx[1][i] * (dLy_dy/wna_)) * tARLQ_ / dens_;
    //     // //}else{
    //     //         // LLx = -(dphi_dx[0][i] * (dLx_dx/data.fWeightFunction(index) - duna_dx*duna_dx - dvna_dx*duna_dy) + 
    //     //         //         dphi_dx[1][i] * (dLx_dy/data.fWeightFunction(index) - duna_dy*duna_dx - dvna_dy*duna_dy)) * tARLQ_ / dens_;
    //     //         // LLy = -(dphi_dx[0][i] * (dLy_dx/data.fWeightFunction(index) - duna_dx*dvna_dx - dvna_dx*dvna_dy) + 
    //     //         //         dphi_dx[1][i] * (dLy_dy/data.fWeightFunction(index) - duna_dy*dvna_dx - dvna_dy*dvna_dy)) * tARLQ_ / dens_;  
    //     //         // LLx = -(dphi_dx[0][i] * (dLx_dx/data.fWeightFunction(index) - dp_dxx) + 
    //     //         //         dphi_dx[1][i] * (dLx_dy/data.fWeightFunction(index) - dp_dxy)) * tARLQ_ / dens_;
    //     //         // LLy = -(dphi_dx[0][i] * (dLy_dx/data.fWeightFunction(index) - dp_dxy) + 
    //     //         //         dphi_dx[1][i] * (dLy_dy/data.fWeightFunction(index) - dp_dyy)) * tARLQ_ / dens_;  

    //     //         // if (iTimeStep > 10){
    //     //         //     LLx +=  -(dphi_dx[0][i]*dax_dx + dphi_dx[1][i]*dax_dy) * tARLQ_;
    //     //         //     LLx +=  -(dphi_dx[0][i]*day_dx + dphi_dx[1][i]*day_dy) * tARLQ_ ;
    //     //         // };

    //     //     //}
            
    //     //     // arlequinStabVector[2*i  ] += (Amx + LLx) * weight_ * djac_ * wna_;
    //     //     // arlequinStabVector[2*i+1] += (Amy + LLy) * weight_ * djac_ * wna_;
    //     // };
    // }
    





};



void CouplingGlobal::ArlequinStabResidual(int &index, MatrixDouble &dphi_dx, VecDouble &phiGlobal, MatrixDouble &dphi_dxGlobal, double &weight_, double &djac_, std::vector<VecDouble> &Rhs){
    
    // auto force = this->Mesh()->getProblemParameters().getForcingFunction();    
    // VecDouble xna_(tshape::Dimension);
    // VecDouble forcingF(tshape::Dimension);
    // auto elglobal = fMeshVector[0]->ElementVec()[fGlobalIndex];
    // xna_ = elglobal->getIntegPointCoordinatesValue(index);
    // // for (int i = 0; i < DIM; i++) xna_[i] = elglobal->setIntegPointWeightFunction(index,i);
    // if (force) force(xna_,forcingF);

    // double WJ = djac_ * weight_ * fMeshVector[0]->ElementVec()[fGlobalIndex]->getIntegPointWeightFunction(index); 

    // //Lagrange Multiplier
    // VecDouble lagM_(this->Mesh()->NState()+1);
    // this->interpolateSolution(index, lagM_);

    // //Lagrange Multiplier Derivatives
    // MatrixDouble dL_dx(this->Mesh()->NState()+1,tshape::Dimension);
    // this->interpolateSolDerivatives(dL_dx);
    // int DIM = tshape::Dimension;

    // if (this->Mesh()->getProblemParameters().ProbType() == EPoisson){
    //     for (int i = 0; i < nphi; i++){
    //         //ARLEQUIN STABILIZATION TERMS
    //         double LLx = 0.;
    //         double LF = 0.;
    //         double shapeFi = data.fPhi[i];
    //         // for (int m = DIM; m--; ) LLx -= dphi_dxGlobal(i,m) * dL_dx(0,m);
    //         // for (int m = DIM; m--; ) LF +=  dphi_dxGlobal(i,m) * forcingF[0] * fMeshVector[0]->ElementVec()[fGlobalIndex]->getIntegPointWeightFunction(index);
    //         // double LF2 = xna_[0]*xna_[0]*elglobal->getIntegPointWeightFunction(index)*shapeFi;

    //         // Rhs[1][i] += (LF2) * weight_ * djac_;
    //         Rhs[1][nphi+i] += (LLx + LF) * weight_ * djac_;
    //     };    
    // } else if (this->Mesh()->getProblemParameters().ProbType() == EElastic) {
    //     VecDouble fieldForce = this->Mesh()->getProblemParameters().GetFieldForce();
    //     MatrixDouble matD(3,2*nphi);
    //     matD.setZero();
    //     MatrixDouble Hooke(3,3);
    //     Hooke.setZero();
    //     // For EPT
    //     double elastic_ = this->Mesh()->getProblemParameters().GetYoungModulus();
    //     double poisson_ = this->Mesh()->getProblemParameters().GetPoissonRatio();
    //     double k = elastic_ / (1. - poisson_ * poisson_);
    //     Hooke(0,0) = k;
    //     Hooke(0,1) = k * poisson_;
    //     Hooke(1,0) = k * poisson_;
    //     Hooke(1,1) = k;
    //     Hooke(2,2) = k * (1. - poisson_) * 0.5;
    //     // Hooke.setIdentity();
    //     for (int j = 0; j < nphi; j++){
    //         matD(0,DIM*j  ) = dphi_dxGlobal(j,0);
    //         matD(1,DIM*j+1) = dphi_dxGlobal(j,1);
    //         matD(2,DIM*j  ) = dphi_dxGlobal(j,1);
    //         matD(2,DIM*j+1) = dphi_dxGlobal(j,0);
    //     }
        
    //     VecDouble strain(3);
    //     strain.setZero();
    //     strain[0] = dL_dx(0,0);
    //     strain[1] = dL_dx(1,1);
    //     strain[2] = dL_dx(0,1)+dL_dx(1,0);
    //     auto aux = matD.transpose() * Hooke * strain * WJ;
    //     // Rhs[1] -= matD.transpose() * Hooke * strain * WJ;

    //     VecDouble forcingF(2);
    //     VecDouble x_ = fMeshVector[0]->ElementVec()[fGlobalIndex]->getIntegPointCoordinatesValue(index);
    //     if (force) force(x_,forcingF);

    //     for (int i = nphi; i--; ){
    //         Rhs[1][2*nphi+2*i  ] += aux[2*i];
    //         Rhs[1][2*nphi+2*i+1] += aux[2*i+1];
    //         double shapeFi = phiGlobal[i];
    //         //External force
    //         double Fx = (fieldForce[0] + forcingF[0]) * shapeFi;
    //         double Fy = (fieldForce[1] + forcingF[1]) * shapeFi;
    //         Rhs[1][2*nphi+2*i  ] += Fx * WJ;
    //         Rhs[1][2*nphi+2*i+1] += Fy * WJ;
    //     };
    // } else {
    //     return;
    //     PanicButton();
    // }
};


