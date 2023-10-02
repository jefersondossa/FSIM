#include "Assemble.h"
#include "ElCouplingLocal.h"
#include "ElCouplingGlobal.h"

void Assemble::Monomodel(Analysis *fAnalysis, int mesh, int64_t startDOF){
    std::cout << "Assembling..." << std::endl;
    int nLocDOF = fAnalysis->MeshVector()[mesh]->NLocDOF(); 
    int nElNodes = fAnalysis->MeshVector()[mesh]->NElNodes(); 
    for (int jel = 0; jel < fAnalysis->MeshVector()[mesh]->NElements(); jel++){   
        if (fAnalysis->MeshVector()[mesh]->part_elem[jel] == 0) {
            //Compute Element matrix
            VecInt connec = fAnalysis->MeshVector()[mesh]->ElementVec()[jel] -> getConnectivity();

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
                            int dof_i = startDOF + nstatei * connec[i] + istate;
                            int dof_j = startDOF + nstatej * connec[j] + jstate; 
                            MatSetValues(fAnalysis->Stiffness(), 1, &dof_i,1, &dof_j, &matrix(nstatei*i+istate,nstatej*j+jstate), ADD_VALUES);
                        }
                    }
                };
                
                //Rhs vector
                for (int istate = 0; istate < nstatei; istate++){
                    int dof_i = startDOF + nstatei * connec[i] + istate;
                    VecSetValues(fAnalysis->Rhs(), 1, &dof_i, &rhs[nstatei*i+istate], ADD_VALUES);
                }
            };
        };
    }; //Elements
    
    

}
void Assemble::Coupling(Analysis *fAnalysis, int64_t startDOF){
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int DIM = fAnalysis->MeshVector()[2]->Dimension();
    int DEG = fAnalysis->MeshVector()[2]->GetDefaultOrder();
    int64_t GloDOF = fAnalysis->MeshVector()[0]->NGlobalDOF();
    int64_t LocDOF = fAnalysis->MeshVector()[1]->NGlobalDOF();

    //Lagrange Multipliers
    for (int jelc=0; jelc< fAnalysis->MeshVector()[2]->NElements(); jelc++){
        auto *elclocal = dynamic_cast<ElCouplingLocal *> (fAnalysis->MeshVector()[2]->ElementVec()[jelc]);
        auto *elcglobal = dynamic_cast<ElCouplingGlobal *> (fAnalysis->MeshVector()[2]->ElementVec()[jelc]);
        int64_t jelcoupled = 0;

        //Determine if the coupling element is from global or local model
        if (!elclocal && !elcglobal) PanicButton();
        if(elclocal){
            jelcoupled = elclocal->GetLocalIndex();
        } else if (elcglobal){
            jelcoupled = elcglobal->GetGlobalIndex();
        } else {
            PanicButton();
        }
        // if (fAnalysis->MeshVector()[1]->part_elem[jel] == rank) {
            
        int nElNodes = fAnalysis->MeshVector()[2]->NElNodes();
        int nLocDOF = fAnalysis->MeshVector()[2]->NLocDOF();

        MatrixDouble matrix(nLocDOF,nLocDOF);
        matrix.setZero();
        VecDouble rhs(2*nLocDOF);
        rhs.setZero();
        auto connecL = fAnalysis->MeshVector()[2]->ElementVec()[jelc]->getConnectivity();
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
        int nState = fAnalysis->MeshVector()[0]->NState(); 
        // if (fAnalysis->MeshVector()[0]->ProbType() == EPoisson) nState = 1;

        //Disperse local contributions into the global matrix
        for (int i = 0; i < nElNodes; i++){
            for (int j = 0; j < nElNodes; j++){
                //COUPLING OPERATOR
                for (int istate = 0; istate < nState; istate++){
                    for (int jstate = 0; jstate < nState; jstate++){
                        if (fabs(matrix(nState*i+istate,nState*j+jstate)) >= 1.e-15){
                            int d_i = GloDOF + LocDOF + nState*connecL[i] + istate;
                            int d_j = startDOF + nState*connec[j] + jstate;
                            MatSetValues(fAnalysis->Stiffness(),1,&d_i,1,&d_j,&matrix(nState*i+istate,nState*j+jstate),ADD_VALUES);
                            MatSetValues(fAnalysis->Stiffness(),1,&d_j,1,&d_i,&matrix(nState*i+istate,nState*j+jstate),ADD_VALUES);
                        };
                    }
                }
                
                // // ARLEQUIN STABILIZATION
                // if (fabs(ArlequinA2(i,j)) >= 1.e-15){
                //     int dof_i = GloDOF + LocDOF + connecL[i];
                //     int dof_j = GloDOF + connec[j];
                //     MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(i,j),ADD_VALUES);
                // };
                // if (fabs(ArlequinA1(i,j)) >= 1.e-15){
                //     int dof_i = GloDOF + LocDOF + connecL[i];
                //     int dof_j = GloDOF + LocDOF + connecL[j];
                //     MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(i,j),ADD_VALUES);
                // };                   
            };
            //RHS VECTOR
            //COUPLING OPERATOR
            for (int istate = 0; istate < nState; istate++){
                int dof_i = startDOF + nState*connec[i]+istate;
                VecSetValues(fAnalysis->Rhs(),1,&dof_i,&rhs[nState*i+istate],ADD_VALUES);
                
                dof_i = GloDOF + LocDOF + nState*connecL[i]+istate;
                VecSetValues(fAnalysis->Rhs(),1,&dof_i,&rhs[nLocDOF+nState*i+istate],ADD_VALUES);
            }
            // dof_i = GloDOF + connec[i];
            // VecSetValues(b,1,&dof_i,&rhsLagMult2[i],ADD_VALUES);

            // //ARLEQUIN STABILIZATION
            // dof_i = GloDOF + LocDOF + connecL[i];
            // VecSetValues(b,1,&dof_i,&RhsArlequin2[i],ADD_VALUES);
        };





            // VecInt connecC;
            // VecInt connec = fMeshVector[1]->ElementVec()[jel] -> getConnectivity();
            // VecInt connecL = glueZoneFine_[l] -> getConnectivity();
            
            // //FINE MESH
            // //Matrices
            // MatrixDouble matC0(dofMatrices,dofMatrices),matC1(dofMatrices,dofMatrices);
            // MatrixDouble matA0(dofMatrices,dofMatrices),matA1(dofMatrices,dofMatrices);
            // MatrixDouble localMV_mat(dofMatrices,dofMatrices);
            // MatrixDouble matE0(dofMatrices,dofMatrices),matE1(dofMatrices,dofMatrices);
            // matC1.setZero();
            // matA1.setZero();
            // localMV_mat.setZero();
            // matE1.setZero();

            // double tARLQ0_, tARLQ1_;
            
            // //Vectors
            // VecDouble vecC0(dofMatrices),vecC1(dofMatrices);
            // VecDouble vecU0(dofMatrices),vecU1(dofMatrices);
            // VecDouble RhsA0(dofMatrices),RhsA1(dofMatrices);
            // VecDouble localMV_vec(dofMatrices);
            // VecDouble vecE0(dofMatrices),vecE1(dofMatrices);
            // vecC1.setZero();
            // vecU1.setZero();
            // localMV_vec.setZero();
            // vecE1.setZero();

            // // FINE MESH
            // //Computes element matrix
            // fMeshVector[1]->ElementVec()[jel] -> getLagrangeMultipliersSameMesh(matC1, vecC1, vecU1);
            
            // if (fArlequinStab != ArlequinStabType::ENoStab){
            //     if (fProbType == ENavierStokes || fProbType == EStokes){
            //         //PSPG and SUPG stabilizations
            //         fMeshVector[1]->ElementVec()[jel] -> getLagrangeMultipliersSUPG_PSPG_SameMesh(localMV_mat,localMV_vec);
            //     }
            //     //Arlequin Stabilization
            //     fMeshVector[1]->ElementVec()[jel] -> getLagrangeMultipliersArlequinSameMesh(matE1, matA1, vecE1);
            // }

            // //COAESE MESH
            // //Counts number of coarse mesh intersecting the fine element
            // int numberIntPoints = fMeshVector[1]->ElementVec()[jel] -> 
            //     getNumberOfIntegrationPoints();
            // int aux;
            
            // std::vector<int> ele, diffElem;
            // ele.clear();
            // diffElem.clear();

            // ele.reserve(3);
            // for (int i=0; i<numberIntPoints; i++){
            //     aux = fMeshVector[1]->ElementVec()[jel] -> 
            //         getIntegPointCorrespondenceElement(i);
            //     ele.push_back(aux);
            //     //std::cout << "Num elem inters " << jel << " " << aux << std::endl;
            // };
            
            // int numElemIntersect = 1;
            // int flag = 0;
            // diffElem.push_back(ele[0]);
            
            // for (int i = 1; i<numberIntPoints; i++){
            //     flag = 0;
            //     for (int j = 0; j<numElemIntersect; j++){
            //         if (ele[i] == diffElem[j]) {
            //             break;
            //         }else{
            //             flag++;
            //         };
            //         if(flag == numElemIntersect){
            //             numElemIntersect++;
            //             diffElem.push_back(ele[i]);
            //         };
            //     };
            // };
            // //Compute the Lagrange Multiplier element matrix
            // for (int ielem = 0; ielem < numElemIntersect; ielem++){
            //     int nElNodes = fMeshVector[1]->NElNodes();
            //     int iElemCoarse = diffElem[ielem];
            //     double pspg = 0;//(*fMeshVector[0]->ElementVec()[iElemCoarse]) -> getPSPG();
            //     VecDouble press_(nElNodes), velX_(nElNodes), velY_(nElNodes), velXPrev_(nElNodes), velYPrev_(nElNodes);

            //     connecC = fMeshVector[0]->ElementVec()[iElemCoarse] -> getConnectivity();

            //     for (int k = 0; k < nElNodes; k++){
            //         press_[k] = fMeshVector[0]->NodeVec()[connecC[k]] -> getPressure();
            //         velX_[k] = fMeshVector[0]->NodeVec()[connecC[k]] -> getVelocity(0);
            //         velY_[k] = fMeshVector[0]->NodeVec()[connecC[k]] -> getVelocity(1);
            //         velXPrev_[k] = fMeshVector[0]->NodeVec()[connecC[k]] -> getPreviousVelocity(0);
            //         velYPrev_[k] = fMeshVector[0]->NodeVec()[connecC[k]] -> getPreviousVelocity(1);
            //     }
                
            //     matC0.setZero();
            //     localMV_mat.setZero();
            //     matE0.setZero();
                
            //     //Vectors
            //     vecC0.setZero();
            //     vecU0.setZero();
            //     localMV_vec.setZero();
            //     vecE0.setZero();
            //     matA0.setZero();
                
            //     fMeshVector[1]->ElementVec()[jel] -> getLagrangeMultipliersDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,velXPrev_,velYPrev_,matC0,vecC0,vecU0);

            //     if (fArlequinStab != ArlequinStabType::ENoStab){
            //         if (fProbType == ENavierStokes || fProbType == EStokes){
            //             fMeshVector[1]->ElementVec()[jel] -> getLagrangeMultipliersSUPG_PSPG_DifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,localMV_mat,localMV_vec);
            //         }

            //         fMeshVector[1]->ElementVec()[jel] -> getLagrangeMultipliersArlequinDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,matE0,matA0,vecE0);
            //     }
            //     MatrixDouble matE = matE1;
            //     stabilizeArlequin(matA0,matA1,matC0,matC1,matE0,vecE0,vecE1,tARLQ0_,tARLQ1_);

            //     matE0 *= tARLQ0_;
            //     matA0 *= tARLQ0_;
            //     vecE0 *= tARLQ0_;
                
            //     switch (fProbType)
            //     {
            //     case ENavierStokes:
            //     case EStokes:
            //         setMatVecValuesLagMultFineCoarse(matC0, localMV_mat, matE0, matA0, 
            //                                      vecU0, vecC0, localMV_vec, vecE0,
            //                                      fMeshVector[0]->ElementVec()[iElemCoarse] -> getConnectivity(), 
            //                                      glueZoneFine_[l] -> getConnectivity());
            //         break;
            //     case EPoisson:
            //         setMatVecValuesLagMultFineCoarsePoisson(matC0, localMV_mat, matE0, matA0, 
            //                                             vecU0, vecC0, localMV_vec, vecE0,
            //                                             fMeshVector[0]->ElementVec()[iElemCoarse] -> getConnectivity(), 
            //                                             glueZoneFine_[l] -> getConnectivity());
            //         break;
            //     case EElastic:
            //         setMatVecValuesLagMultFineCoarseElasticity(matC0, localMV_mat, matE0, matA0, 
            //                                      vecU0, vecC0, localMV_vec, vecE0,
            //                                      fMeshVector[0]->ElementVec()[iElemCoarse] -> getConnectivity(), 
            //                                      glueZoneFine_[l] -> getConnectivity());
            //         break;
                
            //     default:
            //         PanicButton();
            //         break;
            //     }
            // }; //Number of intersections

            // matE1 *= tARLQ1_;
            // matA1 *= tARLQ1_;
            // vecE1 *= tARLQ1_;

            // switch (fProbType)
            // {
            // case ENavierStokes:
            // case EStokes:
            //     setMatVecValuesLagMultFineFine(matC1,localMV_mat,matE1,matA1, 
            //                             vecU1,vecC1,localMV_vec,vecE1,
            //                             fMeshVector[1]->ElementVec()[jel] -> getConnectivity(),
            //                             glueZoneFine_[l] -> getConnectivity());
            //     break;
            // case EPoisson:
            //     setMatVecValuesLagMultFineFinePoisson(matC1,localMV_mat,matE1,matA1, 
            //                                           vecU1,vecC0,localMV_vec,vecE1,
            //                                           fMeshVector[1]->ElementVec()[jel] -> getConnectivity(),
            //                                           glueZoneFine_[l] -> getConnectivity());
            //     break;
            // case EElastic:
            //     setMatVecValuesLagMultFineFineElasticity(matC1,localMV_mat,matE1,matA1, 
            //                                           vecU1,vecC0,localMV_vec,vecE1,
            //                                           fMeshVector[1]->ElementVec()[jel] -> getConnectivity(),
            //                                           glueZoneFine_[l] -> getConnectivity());
            //     break;
            
            // default:
            //     PanicButton();
            //     break;
            // }

            
           
        // }; // if element belongs to the glue zone
    }; // Glue zone
}

void Assemble::Arlequin(Analysis *fAnalysis){
    std::cout << "Assembling Global model..." << std::endl;
    int64_t numDOFGlobal = fAnalysis->MeshVector()[0]->NGlobalDOF();
    Monomodel(fAnalysis,0,0);

    std::cout << "Assembling Local model..." << std::endl;
    int64_t numDOFLocal = fAnalysis->MeshVector()[1]->NGlobalDOF();
    Monomodel(fAnalysis,1,numDOFGlobal);

    // MatView(fAnalysis->Stiffness(),PETSC_VIEWER_STDOUT_WORLD);
    // MatView(fAnalysis->Stiffness(),PETSC_VIEWER_DRAW_WORLD);
    // VecView(fAnalysis->Rhs(),PETSC_VIEWER_STDOUT_WORLD);

    std::cout << "Assembling Coupling operator..." << std::endl;
    Coupling(fAnalysis,numDOFGlobal+numDOFLocal);
}