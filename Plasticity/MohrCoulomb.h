#ifndef MOHRCOULOMB_H
#define MOHRCOULOMB_H

#include "PlasticityModel.h"


class MohrCoulomb : public PlasticityModel
{
private:
    // Internal friction angle
    REAL fInternalFriction;

public:
    /// @brief Linear hardening plasticity model constructor
    /// @param elast elasticity model
    /// @param cohesion cohesion 
    /// @param intfriction Internal friction
    MohrCoulomb(WeakForm *elast, REAL intfriction);

    /// @brief Overloads the updated weak form stiffness matrix computation for the plasticity model
    /// @param index integration point index
    /// @param data integration point data
    /// @param Stiffness vector of stiffness matrices
    void ComputeTangentStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness, Tensor3D &Stress) override;
    
    /// @brief Computes the element error. The exact solution shoul be provided.
    /// @param data integration point data
    /// @param errors vector storing all errors
    void ComputeError(IntPointData &data, VecDouble &errors) override;

    REAL YieldFunction(int &index, IntPointData &data, Tensor3D &Stress) override;

    REAL PlasticMultiplier(int &index, IntPointData &data, Tensor3D &Stress) override;

    void UpdateStateVariables(int &index, IntPointData &data, Tensor3D &Stress) override;
};



#endif