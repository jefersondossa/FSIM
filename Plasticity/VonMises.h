#ifndef VonMises_H
#define VonMises_H

#include "PlasticityModel.h"

/// @brief Implements the linear hardening plasticity model
class VonMises : public PlasticityModel
{
protected:
    double fVonMisesStress;
    double fDevNorm;

    MatrixDouble fMatP;
    Tensor fFlowVector;
    Tensor fDeviatory;

public:
    /// @brief Linear hardening plasticity model constructor
    /// @param elast elasticity model
    /// @param hardModulus hardening modulus
    /// @param yield yield stress
    VonMises(WeakForm *elast);

    /// @brief Overloads the updated weak form stiffness matrix computation for the plasticity model
    /// @param index integration point index
    /// @param data integration point data
    /// @param Stiffness vector of stiffness matrices
    void ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor &Stress) override;
    
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

    double YieldFunction(int &index, IntPointData &data, Tensor &Stress) override;

    Tensor FlowVector(Tensor &Stress);

    double PlasticMultiplier(int &index, IntPointData &data, Tensor &Stress) override;

    void UpdateStateVariables(int &index, IntPointData &data, Tensor &Stress) override;
};



#endif