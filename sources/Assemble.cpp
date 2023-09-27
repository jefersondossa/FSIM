#include "Assemble.h"

// Assemble::Assemble(Analysis *an){
//     fAnalysis = an;
// }

// Assemble::~Assemble()
// {
// }

void Assemble::Monomodel(Analysis *fAnalysis){

    int nLocDOF = fAnalysis->MeshVector()[0]->nLocDOF; 
    int nElNodes = fAnalysis->MeshVector()[0]->nElNodes; 
    for (int jel = 0; jel < fAnalysis->MeshVector()[0]->ElementVec().size(); jel++){   
        if (fAnalysis->MeshVector()[0]->part_elem[jel] == 0) {
            //Compute Element matrix
            VecInt connec = fAnalysis->MeshVector()[0]->ElementVec()[jel] -> getConnectivity();

            MatrixDouble matrix(nLocDOF,nLocDOF);
            matrix.setZero();
            VecDouble rhs(nLocDOF);
            rhs.setZero();

            fAnalysis->MeshVector()[0]->ElementVec()[jel] -> ComputeElContribution(matrix,rhs);

            //Disperse local contributions into the global matrix
            //Stiffness matrix
            for (int i=0; i<nElNodes; i++){
                int nstatei = fAnalysis->MeshVector()[0]->nodes_[connec[i]]->GetNStateVariables();
                for (int j=0; j<nElNodes; j++){
                    int nstatej = fAnalysis->MeshVector()[0]->NodeVec()[connec[j]]->GetNStateVariables();
                    for (int istate = 0; istate < nstatei; istate++){
                        for (int jstate = 0; jstate < nstatej; jstate++){
                            int dof_i = nstatei * connec[i] + istate;
                            int dof_j = nstatej * connec[j] + jstate; 
                            MatSetValues(fAnalysis->Stiffness(), 1, &dof_i,1, &dof_j, &matrix(nstatei*i+istate,nstatej*j+jstate), ADD_VALUES);
                        }
                    }
                };
                
                //Rhs vector
                for (int istate = 0; istate < nstatei; istate++){
                    int dof_i = nstatei * connec[i] + istate;
                    VecSetValues(fAnalysis->Rhs(), 1, &dof_i, &rhs[nstatei*i+istate], ADD_VALUES);
                }
            };
        };
    }; //Elements
    
    PetscErrorCode    ierr;

    //Assemble matrices and vectors
    ierr = MatAssemblyBegin(fAnalysis->Stiffness(),MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(fAnalysis->Stiffness(),MAT_FINAL_ASSEMBLY);
    
    ierr = VecAssemblyBegin(fAnalysis->Rhs());
    ierr = VecAssemblyEnd(fAnalysis->Rhs());

    MatView(fAnalysis->Stiffness(),PETSC_VIEWER_STDOUT_WORLD);
    // MatView(fAnalysis->Stiffness(),PETSC_VIEWER_DRAW_WORLD);
    VecView(fAnalysis->Rhs(),PETSC_VIEWER_STDOUT_WORLD);

}