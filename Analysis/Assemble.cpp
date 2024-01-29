#include "Assemble.h"
#include "CouplingLocal.h"
#include "CouplingGlobal.h"
#include "Arlequin.h"

void Assemble::Monomodel(Analysis *fAnalysis, int mesh, int64_t startDOF){
    std::cout << "Assembling..." << std::endl;

    for (int jel = fAnalysis->MeshVector()[mesh]->NElements(); jel-- ;){   
        if (fAnalysis->MeshVector()[mesh]->part_elem[jel] == 0) {
            //Compute Element matrix
            Element* el = fAnalysis->MeshVector()[mesh]->ElementVec()[jel];
            VecInt connec = el -> getConnectivity();
            int nLocDOF = el->NLocDOF(); 
            int nElNodes = el->NElNodes(); 

            if (nLocDOF == 0) continue;

            MatrixDouble matrix(nLocDOF,nLocDOF);
            matrix.setZero();
            VecDouble rhs(nLocDOF);
            rhs.setZero();

            fAnalysis->MeshVector()[mesh]->ElementVec()[jel] -> ComputeElContribution(matrix,rhs);

            //Disperse local contributions into the global matrix
            //Stiffness matrix
            for (int i=0; i<nElNodes; i++){
                int nstatei = fAnalysis->MeshVector()[mesh]->NodeVec()[connec[i]]->GetNStateVariables();
                for (int j=0; j<nElNodes; j++){
                    int nstatej = fAnalysis->MeshVector()[mesh]->NodeVec()[connec[j]]->GetNStateVariables();
                    for (int istate = 0; istate < nstatei; istate++){
                        for (int jstate = 0; jstate < nstatej; jstate++){
                            int64_t dof_i = startDOF + nstatei * connec[i] + istate;
                            int64_t dof_j = startDOF + nstatej * connec[j] + jstate; 
                            fAnalysis->GlobalMatrix()->AddValueMatrix(dof_i,dof_j,matrix(nstatei*i+istate,nstatej*j+jstate));
                        }
                    }
                };
                
                //Rhs vector
                for (int istate = 0; istate < nstatei; istate++){
                    int64_t dof_i = startDOF + nstatei * connec[i] + istate;
                    fAnalysis->GlobalMatrix()->AddValueRhs(dof_i,rhs[nstatei*i+istate]);
                }
            };
        };
        // std::cout << "Element << " << jel << ", Type = " << fAnalysis->MeshVector()[mesh]->ElementVec()[jel]->PrintType() << std::endl;
        // fAnalysis->GlobalMatrix()->PrintRhs();
    }; //Elements
    
    

}
void Assemble::Coupling(Analysis *fAnalysis, int64_t startDOF){
    
#ifdef HAS_PETSC
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
#endif
    int DIM = fAnalysis->MeshVector()[2]->Dimension();
    int DEG = fAnalysis->MeshVector()[2]->GetDefaultOrder();
    int64_t GloDOF = fAnalysis->MeshVector()[0]->NGlobalDOF();
    int64_t LocDOF = fAnalysis->MeshVector()[1]->NGlobalDOF();

    //Lagrange Multipliers
    for (int64_t jelc=0; jelc< fAnalysis->MeshVector()[2]->NElements(); jelc++){
        auto connecL = fAnalysis->MeshVector()[2]->ElementVec()[jelc]->getConnectivity();
        auto *elclocal = dynamic_cast<CouplingLocal *> (fAnalysis->MeshVector()[2]->ElementVec()[jelc]->GetWeakForm());
        auto *elcglobal = dynamic_cast<CouplingGlobal *> (fAnalysis->MeshVector()[2]->ElementVec()[jelc]->GetWeakForm());
        int64_t jelcoupled = 0;

        //Determine if the coupling element is from global or local model
#ifdef DEBUG_BUILD
        if (!elclocal && !elcglobal) PanicButton();
#endif
        if(elclocal){
            jelcoupled = elclocal->GetLocalIndex();
        } else if (elcglobal){
            jelcoupled = elcglobal->GetGlobalIndex();
        } else {
            PanicButton();
        }
        // if (fAnalysis->MeshVector()[1]->part_elem[jel] == rank) {
            
        int nElNodes = connecL.size();
        int nLocDOF = nElNodes*fAnalysis->MeshVector()[2]->NState();

        std::vector<MatrixDouble> matrix(1);
        matrix[0].resize(nLocDOF,nLocDOF);
        matrix[0].setZero();
        std::vector<VecDouble> rhs(2);
        rhs[0].resize(2*nLocDOF); rhs[1].resize(2*nLocDOF);
        rhs[0].setZero(); rhs[1].setZero();
        VecInt connec;
        int64_t startDOF = 0;
        //Gets the right connectivity
        if(elclocal){
            connec = fAnalysis->MeshVector()[1]->ElementVec()[jelcoupled]->getConnectivity();
            startDOF = GloDOF;
        } else if (elcglobal){
            connec = fAnalysis->MeshVector()[0]->ElementVec()[jelcoupled]->getConnectivity();
        } else {
            PanicButton();
        }
    
        fAnalysis->MeshVector()[2]->ElementVec()[jelc] -> ComputeElContribution(matrix,rhs);
        int nState = fAnalysis->MeshVector()[2]->NState(); 
        int nStateMonomodel = nState;
        // if(fAnalysis->MeshVector()[2]->getProblemParameters().ProbType() == EStokes || fAnalysis->MeshVector()[2]->getProblemParameters().ProbType() == ENavierStokes) nStateMonomodel++;
        if (fAnalysis->ArlequinModel()->getArlequinStabilization() != ENoStab) stabilizeArlequin(fAnalysis,matrix,rhs,jelc);
        
        //Disperse local contributions into the global matrix
        for (int i = 0; i < nElNodes; i++){
            for (int j = 0; j < nElNodes; j++){
                //The COUPLING OPERATOR
                for (int istate = 0; istate < nState; istate++){
                    for (int jstate = 0; jstate < nState; jstate++){
                        if (fabs(matrix[0](nState*i+istate,nState*j+jstate)) >= 1.e-15){
                            int64_t d_i = GloDOF + LocDOF + nState*connecL[i] + istate;
                            int64_t d_j = startDOF + nStateMonomodel*connec[j] + jstate;
                            fAnalysis->GlobalMatrix()->AddValueMatrix(d_i,d_j,matrix[0](nState*i+istate,nState*j+jstate));
                            fAnalysis->GlobalMatrix()->AddValueMatrix(d_j,d_i,matrix[0](nState*i+istate,nState*j+jstate));
                        };
                    }
                }
                //Arlequin Stabilization
                if (fAnalysis->ArlequinModel()->getArlequinStabilization() != ENoStab){
                    for (int istate = 0; istate < nState; istate++){
                        for (int jstate = 0; jstate < nState; jstate++){
                            //Diagonal term
                            if (fabs(matrix[1](nState*i+istate,nState*j+jstate)) >= 1.e-15){
                                int64_t dof_i = GloDOF + LocDOF + nState*connecL[i] + istate;
                                int64_t dof_j = GloDOF + LocDOF + nState*connecL[j] + jstate;
                                fAnalysis->GlobalMatrix()->AddValueMatrix(dof_i,dof_j,matrix[1](nState*i+istate,nState*j+jstate));
                            };
                            //Non Diagonal term
                            for (int iarl = 2; iarl < matrix.size(); iarl++){
                                if (fabs(matrix[iarl](nState*i+istate,nState*j+jstate)) >= 1.e-15){
                                    int64_t dof_i = GloDOF + LocDOF + nState*connecL[i] + istate;
                                    int64_t dof_j = startDOF + nStateMonomodel*connec[j] + jstate;
                                    fAnalysis->GlobalMatrix()->AddValueMatrix(dof_i,dof_j,matrix[iarl](nState*i+istate,nState*j+jstate));
                                };
                            }
                        }
                    }
                }           
            };
            //RHS VECTOR
            //COUPLING OPERATOR
            for (int istate = 0; istate < nState; istate++){
                int64_t dof_i = startDOF + nStateMonomodel*connec[i]+istate;
                fAnalysis->GlobalMatrix()->AddValueRhs(dof_i,rhs[0][nState*i+istate]);
                
                dof_i = GloDOF + LocDOF + nState*connecL[i]+istate;
                fAnalysis->GlobalMatrix()->AddValueRhs(dof_i,rhs[0][nLocDOF+nState*i+istate]);
            }
            //Arlequin Stabilization
            if (fAnalysis->ArlequinModel()->getArlequinStabilization() != ENoStab){
                for (int istate = 0; istate < nState; istate++){
                    int64_t dof_i = startDOF + nStateMonomodel*connec[i]+istate;
                    fAnalysis->GlobalMatrix()->AddValueRhs(dof_i,rhs[1][nState*i+istate]);
                    
                    dof_i = GloDOF + LocDOF + nState*connecL[i]+istate;
                    fAnalysis->GlobalMatrix()->AddValueRhs(dof_i,rhs[1][nLocDOF+nState*i+istate]);
                }
            }
            // dof_i = GloDOF + connec[i];
            // VecSetValues(b,1,&dof_i,&rhsLagMult2[i],ADD_VALUES);

            // //ARLEQUIN STABILIZATION
            // dof_i = GloDOF + LocDOF + connecL[i];
            // VecSetValues(b,1,&dof_i,&RhsArlequin2[i],ADD_VALUES);
        };

    }; // Glue zone
}

void Assemble::Arlequin(Analysis *fAnalysis){
    std::cout << "Assembling Global model..." << std::endl;
    int64_t numDOFGlobal = fAnalysis->MeshVector()[0]->NGlobalDOF();
    Monomodel(fAnalysis,0,0);

    std::cout << "Assembling Local model..." << std::endl;
    int64_t numDOFLocal = fAnalysis->MeshVector()[1]->NGlobalDOF();
    Monomodel(fAnalysis,1,numDOFGlobal);

    // fAnalysis->GlobalMatrix()->PrintMatrix();
    // fAnalysis->GlobalMatrix()->PrintRhs();

    std::cout << "Assembling Coupling operator..." << std::endl;
    Coupling(fAnalysis,numDOFGlobal+numDOFLocal);
}


void Assemble::stabilizeArlequin(Analysis *fAnalysis, std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs, int64_t &element){
    // return;
    //There are in general three main options for taking the norm of a matrix: the L2, L2 and Linfty norms.
    //In eigen they can be simply obtained by:
    // normL2 = mat.norm(); 
    // normL1 = mat.lpNorm<1>(); 
    // normInfty = mat.lpNorm<Infinity>();
    // if (fAnalysis->MeshVector()[0]->getProblemParameters().ProbType() != EPoisson){
    //     // PanicButton();
    // }


    switch (fAnalysis->ArlequinModel()->getArlequinStabilization())
    {
    case ArlequinStabType::ENoStab:
        {
            break;
        }
    
    case ArlequinStabType::EOption1:
        {
            double normC = Stiffness[0].norm();
            double normE = Stiffness[1].norm();
            double normA = Stiffness[2].norm();
            double normB = Rhs[1].norm();
            if (normB == 0) normB = 1.e-10;
            if (normA == 0) normA = 1.e-10;
            if (normE == 0) normE = 1.e-10;
            if (normC == 0) PanicButton();
            double tArlq = 1.e-0;//std::min({normC/normA, normC/normE, normC/normB});
            Stiffness[1] *= tArlq;
            Stiffness[2] *= tArlq;
            Rhs[1] *= tArlq;
            break;
        }

    case ArlequinStabType::EOption2:
        {
            PanicButton();
            // double normC0 = C0.norm();
            // double normC1 = C1.norm();
            // double normA0 = A0.norm();
            // double normA1 = A1.norm();
            // double normB0 = b0.norm();
            // double normB1 = b1.norm();
            // double normE = E.norm();
            // double aux1 = std::min({normC0/normA0, normC0/normE, normC0/normB0, normC1/normA1, normC1/normE, normC1/normB1});
            // double tArlq0 = aux1;
            // double tArlq1 = aux1;
            break;
        }
    
    case ArlequinStabType::EOption3:
        {
            PanicButton();
            // double normC0 = C0.norm();
            // double normC1 = C1.norm();
            // double normA0 = A0.norm();
            // double normA1 = A1.norm();
            // double normE = E.norm();
            // double normB0 = b0.norm();
            // double normB1 = b1.norm();
            
            // tArlq0 = std::min({normC0/normA0, normC0/normE, normC0/normB0});
            // tArlq1 = std::min({normC1/normA1, normC1/normE, normC1/normB1});

            // double normC = std::min(normC0,normC1);
            // tArlq0 *= normC/normC0;
            // tArlq1 *= normC/normC1;
            break;
        }

    case ArlequinStabType::EOption4:
        {
            // ElCouplingGlobal *elglobal = dynamic_cast<ElCouplingGlobal*> (fAnalysis->MeshVector()[2]->ElementVec()[element]);
            // if (elglobal) return;
            PanicButton();
            double normC = Stiffness[0].norm();
            double normE = Stiffness[1].norm();
            double normA = Stiffness[2].norm();
            double normB = Rhs[1].norm();
            if (normB == 0) normB = 1.e-10;
            if (normA == 0) normA = 1.e-10;
            if (normE == 0) normE = 1.e-10;
            if (normC == 0) PanicButton();
            double tArlq = std::min({normC/normA, normC/normE, normC/normB});
            Stiffness[1] *= tArlq;
            Stiffness[2] *= tArlq;
            Rhs[1] *= tArlq;
            break;
        }

    case ArlequinStabType::EOption5:
        {
            // double normC0 = C0.norm();
            // double normC1 = C1.norm();
            // double normA1 = A1.norm();
            // double normE = E.norm();
            // double normB1 = b1.norm();
            // double normC = std::min(normC0,normC1);
            // tArlq0 = 0.;
            // tArlq1 = std::min({normC/normA1, normC/normE, normC/normB1});
            break;
        }
    default:
        PanicButton();
        break;
    }

}