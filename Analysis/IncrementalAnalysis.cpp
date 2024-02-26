#include "IncrementalAnalysis.h"
#include "VTUGenerator.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <ios>

IncrementalAnalysis::IncrementalAnalysis(CompMesh *cmesh, SolverType stype, int nsteps, std::vector<L2Projection *> &bcinc, double tol, int maxIter) 
    : NonLinearAnalysis(cmesh,stype,tol,maxIter), fNSteps(nsteps){
    fIncrementBC = bcinc;
};


IncrementalAnalysis::IncrementalAnalysis(Arlequin *arl, SolverType stype, int nsteps, std::vector<L2Projection *> &bcinc, double tol, int maxIter) 
    : NonLinearAnalysis(arl,stype,tol,maxIter), fNSteps(nsteps){
    fIncrementBC = bcinc;
};

void IncrementalAnalysis::Run(){

    std::vector<VecDouble > increments(fIncrementBC.size());
    //Define the step size for each BC
    for (int i = 0; i < fIncrementBC.size(); i++){
        auto bcval = fIncrementBC[i]->BCValue();
        increments[i] = bcval/fNSteps;
    }
    
    int iStep = 0;
    while (iStep <= fNSteps){
        std::cout << "Computing step... " << iStep << std::endl;
        
        for (int i = 0; i < fIncrementBC.size(); i++){
            fIncrementBC[i]->BCValue() += increments[i] * iStep;
        }

        NonLinearAnalysis::Run();
        iStep++;
    }
}


void IncrementalAnalysis::Run(std::string filename, std::vector<std::string> &scalnames, std::vector<std::string> &vecnames){
    
    std::vector<VecDouble > increments(fIncrementBC.size());
    std::vector<VecDouble> initialbc(fIncrementBC.size());
    //Define the step size for each BC
    if (fIncrementValue.size() == 0){
        for (int i = 0; i < fIncrementBC.size(); i++){
            auto bcval = fIncrementBC[i]->BCValue();
            initialbc[i] = bcval;
            initialbc[i].setZero();
            increments[i] = bcval/fNSteps;
        }
    } else {
        for (int i = 0; i < fIncrementBC.size(); i++){
            initialbc[i] = fIncrementBC[i]->BCValue();
            increments[i] = fIncrementValue; 
        }
    }
    
    int iStep = 0;
    
    while (iStep < fNSteps){
        std::cout << "Computing step... " << iStep << std::endl;
        
        if (fIncrementValue.size() == 0){
            for (int i = 0; i < fIncrementBC.size(); i++){
                fIncrementBC[i]->BCValue() = initialbc[i] + increments[i] * (iStep+1);
            }
        } else {
            for (int i = 0; i < fIncrementBC.size(); i++){
                fIncrementBC[i]->BCValue() = initialbc[i] + increments[i] * (iStep);
            }
        }
        std::ofstream output("plasticity.txt",std::ios::app);
        output << "\n\n------------------STEP = " << iStep << " ------------------" << std::endl;
        NonLinearAnalysis::Run();
        VTUGenerator::PrintResults(MeshVector()[0],filename,scalnames,vecnames,iStep);
        {
            // std::ofstream output("integral.txt",std::ios::app);
            // std::vector<std::string> integrate = {"Solution","DerivativeX","DerivativeY"};
            // std::set<int> matid = {6};
            // std::map<std::string,VecDouble> result;
            // MeshVector()[0]->Integrate(matid,integrate,result);
            // double E = 1.e10;
            // double nu = 0.48;
            // double sigmaY = 848700./sqrt(3.);
            // double aux = E/((1.+nu)*(1.-2.*nu));
            // double valX = result["DerivativeX"][0];
            // double valY = result["DerivativeY"][0];
            // double disp =  fIncrementBC[0]->BCValue()[1];// * 2. * E/(sigmaY*10.);
            // output << disp << " " << valX << " " << -valY << "\n";
        }
        iStep++;
    }
}

void IncrementalAnalysis::UpdateSolution(){
    
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
                double prevsol = this->MeshVector()[imesh]->NodeVec()[i] -> GetSolution(k);
                this->MeshVector()[imesh]->NodeVec()[i] -> SetPreviousSolution(k,prevsol);
                this->MeshVector()[imesh]->NodeVec()[i] -> IncrementSolution(k,val);
            }
            // std::cout << "Solution mesh " << imesh << " " << i <<" " << val <<"\n";
        };
    }
    this->GlobalMatrix()->ClearSolution();
}