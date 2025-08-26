#ifndef VonMises_H
#define VonMises_H

#include "PlasticityModel.h"

/// @brief Implements the linear hardening plasticity model
class VonMises : public PlasticityModel
{
protected:
    double fVonMisesStress;

    MatrixDouble fMatP;
    Tensor3D fFlowVector;
    Tensor3D fDeviatory;

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
    void ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor3D &Stress) override;
    
    /// @brief Computes the element error. The exact solution shoul be provided.
    /// @param data integration point data
    /// @param errors vector storing all errors
    void ComputeError(IntPointData &data, VecDouble &errors) override;

    double YieldFunction(int &index, IntPointData &data, Tensor3D &Stress) override;

    Tensor3D FlowVector(Tensor3D &Stress);

    double PlasticMultiplier(int &index, IntPointData &data, Tensor3D &Stress) override;

    void UpdateStateVariables(int &index, IntPointData &data, Tensor3D &Stress) override;
};



#endif