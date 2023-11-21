#ifndef ELASTICITYPOSITIONAL2D_H
#define ELASTICITYPOSITIONAL2D_H

#include "WeakForm.h"
/// @brief Implements the elasticity 2D problem in the positional formulation;
/// The material is of Saint-Venant-Kirchhoff type
class ElasticityPositional2D : public WeakForm{
protected:
    // True if plane stress, false if plane strain
    bool fPlaneStress = true;
    // Material Young modulus
    double fYoungModulus;
    // Material Poisson's ratio
    double fPoissonRatio;
public:
    //Positional elasticity 2D constructor
    ElasticityPositional2D(int matid, double young, double poisson, bool planes = true);

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