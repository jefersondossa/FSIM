#ifndef DruckerPrager_H
#define DruckerPrager_H

#include "PlasticityModel.h"

/// @brief Implements the linear hardening plasticity model
class DruckerPrager : public PlasticityModel
{
protected:
    REAL fInternalFriction;
    bool fOuterEdges;
    REAL fEta; // n com charme
    REAL fXi; // E com charme
    REAL fDilatancyAngle; 
    REAL fEtaBar;
    REAL fAlpha;
    REAL fBeta;
    bool fApex = false;
    REAL fUpdatedPressure;
    Tensor3D fTrialDevStrain;
    REAL fK;
    REAL fAlpha2;

public:
    /// @brief Linear hardening plasticity model constructor
    /// @param elast elasticity model
    /// @param hardModulus hardening modulus
    /// @param yield yield stress
    /// @param c cohesion
    /// @param phi internal friction
    /// @param psi dilatacy angle
    /// @param oe outer edges
    
    DruckerPrager(WeakForm *elast, REAL phi, REAL psi, bool oe = true);

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