#include "LinearAnalysis.h"
#include "Assemble.h"

void Analysis::Compute(){
    Assemble::Monomodel(this);
} 

void Analysis::UpdateSolution(){
    VecScatter        ctx;
    Vec               All, Allu;
    PetscErrorCode    ierr;

    ierr = VecDuplicate(Rhs(),&All);
    //Gathers the solution vector to the master process
    ierr = VecScatterCreateToAll(Solution(), &ctx, &All);
    ierr = VecScatterBegin(ctx, Solution(), All, INSERT_VALUES, SCATTER_FORWARD);
    ierr = VecScatterEnd(ctx, Solution(), All, INSERT_VALUES, SCATTER_FORWARD);
    ierr = VecScatterDestroy(&ctx);

    ierr = VecScatterCreateToAll(Rhs(), &ctx, &Allu);
    ierr = VecScatterBegin(ctx, Rhs(), Allu, INSERT_VALUES, SCATTER_FORWARD);
    ierr = VecScatterEnd(ctx, Rhs(), Allu, INSERT_VALUES, SCATTER_FORWARD);
    ierr = VecScatterDestroy(&ctx);
    
    //Updates nodal values
    PetscInt Ione = 1;
    PetscInt Ii;
    PetscScalar val;
        
    for (int i = 0; i < this->Mesh()->NodeVec().size(); ++i){
        int nstate = this->Mesh()->NodeVec()[i]->GetNStateVariables();
        for (int k = 0; k<nstate; k++){
            Ii = nstate*i+k;
            ierr = VecGetValues(All, Ione, &Ii, &val);
            this->Mesh()->NodeVec()[i] -> SetSolution(k,val);
        }
    };
    ierr = VecDestroy(&All); 
    ierr = VecDestroy(&Allu); 
}

