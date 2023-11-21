#ifndef ELEMENTNAVIERSTOKES_H
#define ELEMENTNAVIERSTOKES_H

#include "Stokes.h"

/// @brief Implements the stabilized Navier-Stokes weak form: SUPG, PSPG and LSIC stabilizations are employed
class NavierStokes : public Stokes{
public:
    /// @brief Navier-Stokes flow weak form constructor
    /// @param matid physical tag
    /// @param dim problem dimension (2 or 3)
    /// @param density fluid density
    /// @param viscosity fluid dynamic viscosity
    NavierStokes(int matid, int dim, double density, double viscosity);

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

};


#endif