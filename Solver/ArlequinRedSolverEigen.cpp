#include "ArlequinRedSolverEigen.h"
#include "ArlequinMatRedEigen.h"
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/OrderingMethods>


ArlequinRedSolverEigen::ArlequinRedSolverEigen(Analysis *an) : LinearSolver(an){
    
}

ArlequinRedSolverEigen::~ArlequinRedSolverEigen(){

}

void ArlequinRedSolverEigen::Solve(){
    
    auto * emat = dynamic_cast<ArlequinMatRedEigen*> (fAnalysis->GlobalMatrix());
#ifdef DEBUG_BUILD
    if (!emat) PanicButton();
#endif

    switch (fAnalysis->SType())
    {
    case SolverType::ELU:
        {
            SparseLU<SparseMat,COLAMDOrdering<int>>   solverLambda, solverU0, solverU1;
            solverU0.analyzePattern(emat->K0());
            solverU0.factorize(emat->K0());

            solverU1.analyzePattern(emat->K1());
            solverU1.factorize(emat->K1());
            //Reduced matrix;
            SparseMat aux = emat->E() - emat->L0T()*solverU0.solve(emat->L0()) - emat->L1T()*solverU1.solve(emat->L1());
            //Only needed to avoid singular matrix- It happens when dirichlet BC is applied to nodes in the gluing zone
            for (auto i = 0; i < aux.rows(); i++){
                if (fabs(aux.coeffRef(i,i)) < 1.e-10){
                    aux.coeffRef(i,i) = 1.;
                }
            }
            solverLambda.analyzePattern(aux);
            solverLambda.factorize(aux);
            emat->Lambda() = solverLambda.solve(emat->G()-emat->L0T()*solverU0.solve(emat->F0())-emat->L1T()*solverU1.solve(emat->F1()));
            emat->U0() = solverU0.solve(emat->F0()-emat->L0()*emat->Lambda());
            emat->U1() = solverU1.solve(emat->F1()-emat->L1()*emat->Lambda());
        }
        break;
    case SolverType::ELLt:
        {
            SimplicialLLT<SparseMat>   solverLambda, solverU0, solverU1;
            solverU0.analyzePattern(emat->K0());
            solverU0.factorize(emat->K0());

            solverU1.analyzePattern(emat->K1());
            solverU1.factorize(emat->K1());
            //Reduced matrix;
            SparseMat aux = emat->E() - emat->L0T()*solverU0.solve(emat->L0()) - emat->L1T()*solverU1.solve(emat->L1());
            //Only needed to avoid singular matrix- It happens when dirichlet BC is applied to nodes in the gluing zone
            for (auto i = 0; i < aux.rows(); i++){
                if (fabs(aux.coeffRef(i,i)) < 1.e-10){
                    aux.coeffRef(i,i) = 1.;
                }
            }
            solverLambda.analyzePattern(aux);
            solverLambda.factorize(aux);
            emat->Lambda() = solverLambda.solve(emat->G()-emat->L0T()*solverU0.solve(emat->F0())-emat->L1T()*solverU1.solve(emat->F1()));
            emat->U0() = solverU0.solve(emat->F0()-emat->L0()*emat->Lambda());
            emat->U1() = solverU1.solve(emat->F1()-emat->L1()*emat->Lambda());
        }
        break;
    case SolverType::ELDLt:
        {
            SimplicialLDLT<SparseMat>   solverLambda, solverU0, solverU1;
            solverU0.analyzePattern(emat->K0());
            solverU0.factorize(emat->K0());

            solverU1.analyzePattern(emat->K1());
            solverU1.factorize(emat->K1());
            //Reduced matrix;
            SparseMat aux = emat->E() - emat->L0T()*solverU0.solve(emat->L0()) - emat->L1T()*solverU1.solve(emat->L1());
            //Only needed to avoid singular matrix- It happens when dirichlet BC is applied to nodes in the gluing zone
            for (auto i = 0; i < aux.rows(); i++){
                if (fabs(aux.coeffRef(i,i)) < 1.e-10){
                    aux.coeffRef(i,i) = 1.;
                }
            }
            solverLambda.analyzePattern(aux);
            solverLambda.factorize(aux);
            emat->Lambda() = solverLambda.solve(emat->G()-emat->L0T()*solverU0.solve(emat->F0())-emat->L1T()*solverU1.solve(emat->F1()));
            emat->U0() = solverU0.solve(emat->F0()-emat->L0()*emat->Lambda());
            emat->U1() = solverU1.solve(emat->F1()-emat->L1()*emat->Lambda()); 
        }
        break;
    case SolverType::EQR:
        {
            SparseQR<SparseMat,COLAMDOrdering<int>>   solverLambda, solverU0, solverU1;
            solverU0.analyzePattern(emat->K0());
            solverU0.factorize(emat->K0());

            solverU1.analyzePattern(emat->K1());
            solverU1.factorize(emat->K1());
            //Reduced matrix;
            SparseMat aux = emat->E() - emat->L0T()*solverU0.solve(emat->L0()) - emat->L1T()*solverU1.solve(emat->L1());
            //Only needed to avoid singular matrix- It happens when dirichlet BC is applied to nodes in the gluing zone
            for (auto i = 0; i < aux.rows(); i++){
                if (fabs(aux.coeffRef(i,i)) < 1.e-10){
                    aux.coeffRef(i,i) = 1.;
                }
            }
            solverLambda.analyzePattern(aux);
            solverLambda.factorize(aux);
            emat->Lambda() = solverLambda.solve(emat->G()-emat->L0T()*solverU0.solve(emat->F0())-emat->L1T()*solverU1.solve(emat->F1()));
            emat->U0() = solverU0.solve(emat->F0()-emat->L0()*emat->Lambda());
            emat->U1() = solverU1.solve(emat->F1()-emat->L1()*emat->Lambda());
        }
        break;
    case SolverType::ECG:
        {
            ConjugateGradient<SparseMat>   solverLambda, solverU0, solverU1;
            solverU0.analyzePattern(emat->K0());
            solverU0.factorize(emat->K0());

            solverU1.analyzePattern(emat->K1());
            solverU1.factorize(emat->K1());
            //Reduced matrix;
            SparseMat aux = emat->E() - emat->L0T()*solverU0.solve(emat->L0()) - emat->L1T()*solverU1.solve(emat->L1());
            //Only needed to avoid singular matrix- It happens when dirichlet BC is applied to nodes in the gluing zone
            for (auto i = 0; i < aux.rows(); i++){
                if (fabs(aux.coeffRef(i,i)) < 1.e-10){
                    aux.coeffRef(i,i) = 1.;
                }
            }
            solverLambda.analyzePattern(aux);
            solverLambda.factorize(aux);
            emat->Lambda() = solverLambda.solve(emat->G()-emat->L0T()*solverU0.solve(emat->F0())-emat->L1T()*solverU1.solve(emat->F1()));
            emat->U0() = solverU0.solve(emat->F0()-emat->L0()*emat->Lambda());
            emat->U1() = solverU1.solve(emat->F1()-emat->L1()*emat->Lambda());
        }
        break;
    case SolverType::EBiCGStab:
        {
            BiCGSTAB<SparseMat>   solverLambda, solverU0, solverU1;
            solverU0.analyzePattern(emat->K0());
            solverU0.factorize(emat->K0());

            solverU1.analyzePattern(emat->K1());
            solverU1.factorize(emat->K1());
            //Reduced matrix;
            SparseMat aux = emat->E() - emat->L0T()*solverU0.solve(emat->L0()) - emat->L1T()*solverU1.solve(emat->L1());
            //Only needed to avoid singular matrix- It happens when dirichlet BC is applied to nodes in the gluing zone
            for (auto i = 0; i < aux.rows(); i++){
                if (fabs(aux.coeffRef(i,i)) < 1.e-10){
                    aux.coeffRef(i,i) = 1.;
                }
            }
            solverLambda.analyzePattern(aux);
            solverLambda.factorize(aux);
            emat->Lambda() = solverLambda.solve(emat->G()-emat->L0T()*solverU0.solve(emat->F0())-emat->L1T()*solverU1.solve(emat->F1()));
            emat->U0() = solverU0.solve(emat->F0()-emat->L0()*emat->Lambda());
            emat->U1() = solverU1.solve(emat->F1()-emat->L1()*emat->Lambda());
        }
        break;
    case SolverType::ELSCG:
        {
            LeastSquaresConjugateGradient<SparseMat>   solverLambda, solverU0, solverU1;
            solverU0.analyzePattern(emat->K0());
            solverU0.factorize(emat->K0());

            solverU1.analyzePattern(emat->K1());
            solverU1.factorize(emat->K1());
            //Reduced matrix;
            SparseMat aux = emat->E() - emat->L0T()*solverU0.solve(emat->L0()) - emat->L1T()*solverU1.solve(emat->L1());
            //Only needed to avoid singular matrix- It happens when dirichlet BC is applied to nodes in the gluing zone
            for (auto i = 0; i < aux.rows(); i++){
                if (fabs(aux.coeffRef(i,i)) < 1.e-10){
                    aux.coeffRef(i,i) = 1.;
                }
            }
            solverLambda.analyzePattern(aux);
            solverLambda.factorize(aux);
            emat->Lambda() = solverLambda.solve(emat->G()-emat->L0T()*solverU0.solve(emat->F0())-emat->L1T()*solverU1.solve(emat->F1()));
            emat->U0() = solverU0.solve(emat->F0()-emat->L0()*emat->Lambda());
            emat->U1() = solverU1.solve(emat->F1()-emat->L1()*emat->Lambda());
        }
        break;
    
    default:
        std::cout << "Please select an Eigen supported solver\n";
        PanicButton();
        break;
    }

    
}