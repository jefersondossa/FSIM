#include "NonLinearAnalysis.h"
#include "Assemble.h"
#include "L2Projection.h"

void NonLinearAnalysis::UpdateSolution(){
    
    this->GlobalMatrix()->ExpandSolution();

    //Updates nodal values
    int64_t Ione = 1;
    int64_t Ii;
    double val;
    
    int64_t nstartDOF = 0;
    for (int imesh = 0; imesh < this->MeshVector().size(); imesh++){
        if (imesh > 0) nstartDOF += this->MeshVector()[imesh-1]->NGlobalDOF();
        for (int i = 0; i < this->MeshVector()[imesh]->NNodes(); ++i){
            int nstate = this->MeshVector()[imesh]->NState();
            for (int k = 0; k<nstate; k++){
                Ii = nstartDOF + nstate*i+k;
                val = this->GlobalMatrix()->GetValueSolution(Ii);
                // ierr = VecGetValues(All, Ione, &Ii, &val);
                this->MeshVector()[imesh]->NodeVec()[i] -> IncrementSolution(k,val);
            }
            // std::cout << "Solution mesh " << imesh << " " << i <<" " << val <<"\n";
        };
    }
    this->GlobalMatrix()->ClearSolution();
}

void NonLinearAnalysis::Run(){
    
    double NRL2norm = 1000.;
    int iteration = 0;
    
    //Save Previous solution
    for (int imesh = 0; imesh < this->MeshVector().size(); imesh++){
        for (int64_t inode = 0; inode < this->MeshVector()[imesh]->NNodes(); inode++){
            int nstate = this->MeshVector()[imesh]->NState();
            for (int istate = 0; istate < nstate; istate++){
                double Sol = this->MeshVector()[imesh]->NodeVec()[inode]->GetSolution(istate);
                this->MeshVector()[imesh]->NodeVec()[inode]->SetPreviousSolution(istate,Sol);
            }
        }
    }
    //Iterative Process
    while (NRL2norm > fTolerance && iteration < fMaxIterations)
    {   
        std::clock_t t3 = std::clock();
        Compute();
        std::clock_t t4 = std::clock();
        std::cout << "Time assembling = " << 1000.*(t4-t3)/CLOCKS_PER_SEC/1000. << "s \n";
        Solve();
        std::clock_t t5 = std::clock();
        std::cout << "Time Solving = " << 1000.*(t5-t4)/CLOCKS_PER_SEC/1000. << "s \n";
        UpdateSolution();
        // this->GlobalMatrix()->PrintMatrix();
        // this->GlobalMatrix()->PrintRhs();
        // this->GlobalMatrix()->PrintSolution();
        NRL2norm = this->GlobalMatrix()->SolutionNorm();
        std::cout << "Iteration " << iteration++ << ", Newton-Raphson residual = " << NRL2norm << std::endl;
        this->GlobalMatrix()->ZeroMatrix();
        this->GlobalMatrix()->ZeroRhs();
        this->GlobalMatrix()->ZeroSolution();
        for (int64_t i=0; i<NEquations(); i++){
            double val = 1.e-20;
            this->GlobalMatrix()->AddValueMatrix(i,i,val);
        }
    }
    
    


}
