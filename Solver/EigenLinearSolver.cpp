#include "EigenLinearSolver.h"
#include "EigenSpMatrix.h"
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/OrderingMethods>


EigenLinearSolver::EigenLinearSolver(Analysis *an) : LinearSolver(an){
    
}

EigenLinearSolver::~EigenLinearSolver(){

}

void EigenLinearSolver::Solve(){
    
    auto * emat = dynamic_cast<EigenSpMatrix*> (fAnalysis->GlobalMatrix());
#ifdef DEBUG_BUILD
    if (!emat) PanicButton();
#endif

    switch (fAnalysis->SType())
    {
    case ELU:
        {
            SparseLU<SparseMat,COLAMDOrdering<int>>   solver;
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix());
            emat->Solution() = solver.solve(emat->Rhs()); 
            std::cout << "Mat determinant = " << solver.determinant() << std::endl;
        }
        break;
    case ELLt:
        {
            SimplicialLLT<SparseMat>   solver;
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix());
            emat->Solution() = solver.solve(emat->Rhs()); 
        }
        break;
    case ELDLt:
        {
            SimplicialLDLT<SparseMat>   solver;
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix());
            emat->Solution() = solver.solve(emat->Rhs()); 
        }
        break;
    case EQR:
        {
            SparseQR<SparseMat,COLAMDOrdering<int>>   solver;
            // SparseQR<SparseMat,AMDOrdering<int>>   solver;
            // SparseQR<SparseMat,NaturalOrdering<int>>   solver;
            emat->Matrix().makeCompressed();
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix());
            emat->Solution() = solver.solve(emat->Rhs()); 
        }
        break;
    case ECG:
        {
            ConjugateGradient<SparseMat>   solver;
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix());
            emat->Solution() = solver.solve(emat->Rhs()); 
        }
        break;
    case EBiCGStab:
        {
            BiCGSTAB<SparseMat>   solver;
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix()); 
            emat->Solution() = solver.solve(emat->Rhs()); 
        }
        break;
    case ELSCG:
        {
            LeastSquaresConjugateGradient<SparseMat>   solver;
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix());
            emat->Solution() = solver.solve(emat->Rhs()); 
        }
        break;
    
    default:
        std::cout << "Please select an Eigen supported solver\n";
        PanicButton();
        break;
    }

    
}