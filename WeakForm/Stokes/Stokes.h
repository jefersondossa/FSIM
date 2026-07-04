#ifndef ELEMENTSTOKES_H
#define ELEMENTSTOKES_H

#include "WeakForm.h"

/// @brief Implements the (PSPG) stabilized Stokes equation weak form
class Stokes : public WeakForm{
protected:
    // Fluid dynamic viscosity
    REAL fViscosity;
    // Fluid density
    REAL fDensity;
    // Stabilization parameters
    REAL tPSPG_, tSUPG_, tLSIC_; 

    /// @brief Computes the PSPG, SUPG and LSIC stabilization parameters
    /// @param index integration point index
    /// @param data integration point data
    virtual void GetStabilizationParameter(int &index, IntPointData &data);

public:
    /// @brief Stokes weak form constructor
    /// @param matid the physical tag
    /// @param dim problem dimension (2 or 3)
    /// @param density fluid density
    /// @param viscosity fluid dynamic viscosity
    Stokes(int matid, int dim, REAL density, REAL viscosity);

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
};


#endif