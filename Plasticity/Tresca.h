#ifndef Tresca_H
#define Tresca_H

#include "PlasticityModel.h"
#include "PlasticStep.h"

/// @brief Implements the linear hardening plasticity model
class Tresca : public PlasticityModel
{
private:
    enum ETrescaDirection{EMainPlane, ERightCorner, ELeftCorner};
    //The principal stresses
    double fS1, fS2, fS3;
    MatrixDouble fPrincipalDirections;
    ETrescaDirection fReturnDirection;

public:
    /// @brief Linear hardening plasticity model constructor
    /// @param elast elasticity model
    /// @param hardModulus hardening modulus
    /// @param yield yield stress
    Tresca(WeakForm *elast);

    /// @brief Overloads the updated weak form stiffness matrix computation for the plasticity model
    /// @param index integration point index
    /// @param data integration point data
    /// @param Stiffness vector of stiffness matrices
    void ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor &Stress) override;
    
    double YieldFunction(int &index, IntPointData &data, Tensor &Stress) override;

    Tensor FlowVector(Tensor &Stress);

    double PlasticMultiplier(int &index, IntPointData &data, Tensor &Stress) override;

    void UpdateStateVariables(int &index, IntPointData &data, Tensor &Stress) override;
    
};



#endif