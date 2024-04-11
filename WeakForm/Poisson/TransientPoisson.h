#ifndef TRANSIENTPOISSON_H
#define TRANSIENTPOISSON_H

#include "Poisson.h"
#include "TransientWeakForm.h"

/// @brief Implements the Poisson's equation weak form
class TransientPoisson : public Poisson, public TransientWeakForm{
protected:
    // Damping
    double fDamping;
    // Density
    double fDensity;
    // Integration parameters
    double fBeta = 0.25;
    double fGamma = 0.5;
    double fSpectralRadius = 1.;

public:
    /// @brief Poisson weak for constructor
    /// @param matid physical tag
    /// @param dim problem dimension
    /// @param nState number of state variables
    TransientPoisson(int matid, int dim, double damp, double dens, double dt, TimeIntegScheme integscheme = ENewmark, int nState = 1);

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

    void SetNewmarkParameters(double beta, double gamma){
        fBeta = beta;
        fGamma = gamma;
    }

    void UpdateTimeDerivatives(CompMesh *cmesh) override;

};


#endif