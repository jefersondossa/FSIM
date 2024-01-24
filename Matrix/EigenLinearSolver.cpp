#include "EigenLinearSolver.h"
#include "EigenSpMatrix.h"
#include <Eigen/SparseLU>

EigenLinearSolver::EigenLinearSolver(Analysis *an) : LinearSolver(an){
    
}

EigenLinearSolver::~EigenLinearSolver(){

}

void EigenLinearSolver::Solve(){
    SparseLU<SparseMat, COLAMDOrdering<int> >   solver;
    
    auto * emat = dynamic_cast<EigenSpMatrix*> (fAnalysis->GlobalMatrix());
    if (!emat) PanicButton();

    solver.analyzePattern(emat->Matrix()); 
    // Compute the numerical factorization 
    solver.factorize(emat->Matrix()); 
    //Use the factors to solve the linear system 
    emat->Solution() = solver.solve(emat->Rhs()); 
}