#ifndef MOHRCOULOMB_H
#define MOHRCOULOMB_H

#include "PlasticityModel.h"


class MohrCoulomb : public PlasticityModel
{
private:
    // The updated constitutive matrix
    MatrixDouble fConstitutiveMatrix;

    // Internal friction angle
    double fInternalFriction;

    // Cohesion modulus
    double fCohesion;

    // Young modulus
    double fYoungModulus;

    // Plastic strain tensor
    MatrixDouble fPlasticStrain;

    // Total strain tensor
    MatrixDouble fTotalStrain;

public:
    /// @brief Linear hardening plasticity model constructor
    /// @param elast elasticity model
    /// @param cohesion cohesion 
    /// @param intfriction Internal friction
    MohrCoulomb(WeakForm *elast, double cohesion, double intfriction);

    /// @brief Overloads the updated weak form stiffness matrix computation for the plasticity model
    /// @param index integration point index
    /// @param data integration point data
    /// @param Stiffness vector of stiffness matrices
    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    /// @brief Returns the updated weak form residual vector for the plasticity model 
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