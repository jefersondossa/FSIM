#include "EigenLinearSolver.h"
#include "EigenSpMatrix.h"
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/OrderingMethods>

#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

namespace
{
    void NormalizeSystem(SparseMat &kMatrix, MatrixDouble &rhs)
    {
        for (int k = 0; k < kMatrix.outerSize(); ++k)
        {
            double max_val = 0.0;
            for (SparseMatrix<double>::InnerIterator it(kMatrix, k); it; ++it)
            {
                max_val = std::max(std::fabs(max_val), std::fabs(it.value()));
            }

            for (SparseMatrix<double>::InnerIterator it(kMatrix, k); it; ++it)
            {
                it.valueRef() /= max_val;

                if (std::fabs(it.valueRef()) < 1e-12)
                {
                    it.valueRef() = 0.0;
                }
            }

            rhs.data()[k] /= max_val;

            if (std::fabs(rhs.data()[k]) < 1e-12)
            {
                rhs.data()[k] = 0.0;
            }
        }
    }

    void EigenToCRS(const SparseMat &kMatrix, const MatrixDouble &rhs,
                    std::vector<int> &ptr_k_f, std::vector<int> &col_k_f, std::vector<double> &val_k_f,
                    std::vector<double> &rhs_v)
    {
        ptr_k_f.clear();
        col_k_f.clear();
        val_k_f.clear();
        rhs_v.clear();

        assert(kMatrix.rows() == kMatrix.cols());
        ptr_k_f.push_back(0);
        const auto ndofs = kMatrix.cols();
        ptr_k_f.reserve(ndofs + 1);
        col_k_f.reserve(ndofs * 128);
        val_k_f.reserve(ndofs * 128);

        for (int k = 0; k < kMatrix.outerSize(); ++k)
        {
            for (SparseMatrix<double>::InnerIterator it(kMatrix, k); it; ++it)
            {
                if (it.value() != 0.0)
                {
                    col_k_f.push_back(it.index());
                    val_k_f.push_back(it.value());
                }
            }
            ptr_k_f.push_back(col_k_f.size());
        }

        rhs_v.reserve(rhs.rows());

        for (size_t i = 0; i < rhs.rows(); i++)
        {
            rhs_v.push_back(rhs.data()[i]);
        }
    }

    MatrixDouble CRSSolutionToEigen(const std::vector<double> &v_U_f)
    {
        MatrixDouble solution;
        solution.resize(v_U_f.size(), 1);

        for (std::size_t i = 0; i < v_U_f.size(); i++)
        {
            solution(i) = v_U_f[i];
        }

        return solution;
    }

    MatrixDouble SolveUsingAMGCL(SparseMat &kMatrix, MatrixDouble &rhs)
    {
        using Backend = amgcl::backend::builtin<double>;

        using Preconditioner = amgcl::amg<
            Backend,
            amgcl::coarsening::smoothed_aggregation,
            amgcl::relaxation::spai0>;

        using LSSolver = amgcl::solver::bicgstab<Backend>;

        using Solver = amgcl::make_solver<Preconditioner, LSSolver>;

        constexpr double kTolerance = 1e-6;
        constexpr std::size_t kMaxIterations = 1'000;

        Solver::params prm;
        prm.solver.tol = kTolerance;
        prm.solver.maxiter = kMaxIterations;

        NormalizeSystem(kMatrix, rhs);

        std::vector<int> ptr_k_f;
        std::vector<int> col_k_f;
        std::vector<double> val_k_f;
        std::vector<double> rhs_v;

        EigenToCRS(kMatrix, rhs, ptr_k_f, col_k_f, val_k_f, rhs_v);

#if 0
        const auto get_norm_sparse_infinity = [](const auto& matrix) {
            double k_norm_infinity = std::numeric_limits<double>::min();

            for (int k=0; k<matrix.outerSize(); ++k)
            {
                double curr_norm = 0;
                for (SparseMatrix<double>::InnerIterator it(matrix,k); it; ++it)
                {
                    curr_norm += std::fabs(it.value());
                }

                k_norm_infinity = std::max(curr_norm, k_norm_infinity);
            }

            return k_norm_infinity;
        };

        const auto norm_k_matrix = get_norm_sparse_infinity(kMatrix);

        const auto get_norm_dense_infinity = [](const auto& matrix) {
            double k_norm_infinity = std::numeric_limits<double>::min();

            for (int k=0; k<matrix.rows(); ++k)
            {
                double curr_norm = 0;
                for (int j = 0; j < matrix.cols(); j++)
                {
                    curr_norm += std::fabs(matrix(k,j));
                }

                k_norm_infinity = std::max(curr_norm, k_norm_infinity);
            }

            return k_norm_infinity;
        };

        const auto norm_inverse = get_norm_dense_infinity(MatrixXd(kMatrix).inverse());

        std::cout << "Condition number: " << norm_k_matrix * norm_inverse << std::endl;
#endif

        static std::vector<double> v_U_f{};
        const auto ndofs = kMatrix.cols();
        if(v_U_f.size() != ndofs)
        {
            v_U_f.resize(ndofs);
        }

        Solver solve{std::tie(ndofs, ptr_k_f, col_k_f, val_k_f), prm};

        int nItOut;
        double errorOut;
        std::tie(nItOut, errorOut) = solve(rhs_v, std::move(v_U_f));

        std::cout << "Solved using AMGCL." << std::endl;
        std::cout << "    # iterations: " << nItOut << std::endl;
        std::cout << "    Final error: " << errorOut << std::endl;

        return CRSSolutionToEigen(v_U_f);
    }
}

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
    case SolverType::ELU:
        {
            SparseLU<SparseMat,COLAMDOrdering<int>>   solver;
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix());
            emat->Solution() = solver.solve(emat->Rhs()); 
            std::cout << "Mat determinant = " << solver.determinant() << std::endl;
        }
        break;
    case SolverType::ELLt:
        {
            SimplicialLLT<SparseMat>   solver;
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix());
            emat->Solution() = solver.solve(emat->Rhs()); 
        }
        break;
    case SolverType::ELDLt:
        {
            SimplicialLDLT<SparseMat>   solver;
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix());
            emat->Solution() = solver.solve(emat->Rhs()); 
        }
        break;
    case SolverType::EQR:
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
    case SolverType::ECG:
        {
            ConjugateGradient<SparseMat>   solver;
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix());
            emat->Solution() = solver.solve(emat->Rhs()); 
        }
        break;
    case SolverType::EBiCGStab:
        {
            BiCGSTAB<SparseMat>   solver;
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix()); 
            emat->Solution() = solver.solve(emat->Rhs()); 
        }
        break;
    case SolverType::ELSCG:
        {
            LeastSquaresConjugateGradient<SparseMat>   solver;
            solver.analyzePattern(emat->Matrix());
            solver.factorize(emat->Matrix());
            emat->Solution() = solver.solve(emat->Rhs()); 
        }
        break;
        case SolverType::AMGCLBiCGStab:
        {
            emat->Solution() = SolveUsingAMGCL(emat->Matrix(), emat->Rhs());
        }
        break;

        default:
            std::cout << "Please select an Eigen supported solver\n";
            PanicButton();
            break;
    }

    
}