#ifndef TRANSIENTELASTICITY2D_H
#define TRANSIENTELASTICITY2D_H

#include "Elasticity2D.h"

/// @brief Implements the linear (Hooke's law) transient 2D elasticity problem with newmark integration
class TransientElasticity2D : public Elasticity2D{
public:  
    enum TimeIntegScheme {ENewmark, EGeneralizedAlpha};

protected:
    // Damping
    double fDamping;
    // Density
    double fDensity;
    // Time Step
    double fTimeStep;
    // Integration parameters
    double fBeta = 0.25;
    double fGamma = 0.5;
    double fSpectralRadius = 1.;

    TimeIntegScheme fIntegScheme;

public:  
    // Elasticity 2D class constructor
    TransientElasticity2D(int matid, double young, double poisson, bool planes, double damp, double dens, double dt, TimeIntegScheme integscheme = ENewmark);

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
};


#endif