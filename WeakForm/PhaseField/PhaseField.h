#ifndef POISSON_H
#define POISSON_H

#include "TransientWeakForm.h"
#include "WeakForm.h"

static constexpr double KK = 1.8;

/// @brief Implements the PhaseField's equation weak form
class PhaseField : public TransientWeakForm, public WeakForm{
private:
    static constexpr double dx = 0.01;
    static constexpr double fM = 1.0;
    static constexpr double fGamma = 8*dx * KK;
    static constexpr double fKsi = 1e-4;
    static constexpr double fBeta = 1.5 * KK;
    static constexpr double fDt = 0.25 * (dx * dx) / (4 * fGamma * fKsi) / (6*8*4);//5e-3;

private:
    void ComputeStiffnessStatic(int &index, IntPointData &data, MatrixDouble &Stiffness);

    void ComputeResidualStatic(int &index, IntPointData &data, VecDouble &Rhs);
public:
    /// @brief PhaseField weak for constructor
    /// @param matid physical tag
    /// @param dim problem dimension
    /// @param nState number of state variables
    PhaseField(int matid, int dim);

    /// @brief Overloads the weak form stiffness matrix computation in the case more than one contribution is provided
    /// @param index integration point index
    /// @param data integration point data
    /// @param Stiffness vector of stiffness matrices
    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    /// @brief Returns the weak form residual vector. It should never be called in this class
    /// @param index integration point index
    /// @param data integration point data
    /// @param Rhs residual vector
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override;
    
    /// @brief Computes the element error. The exact solution shoul be provided.
    /// @param data integration point data
    /// @param errors vector storing all errors
    void ComputeError(IntPointData &data, VecDouble &errors) override;

    /// @brief Returns the variable index of a given solution variable
    /// @param name solution variable name
    /// @return solution variable's index
    int VariableIndex(const std::string &name) const override;

    /// @brief Returns the number of solution variables of a given post processing variable
    /// @param var solution variable's index
    /// @return number of solution variables
    int NSolutionVariables(int var) const override;
    
    /// @brief Post process the results for a given solution variable. It should never be called here, but in the derived weak form.
    /// @param data integration point data
    /// @param var solution variable's index
    /// @param Sol solution vector
    void Solution(IntPointData &data, int var, VecDouble &Sol) override;

    void UpdateTimeDerivatives(CompMesh *cmesh) override;
};


#endif