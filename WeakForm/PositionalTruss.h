#ifndef PositionalTruss_H
#define PositionalTruss_H

#include "WeakForm.h"
/// @brief Implements the truss element in positional formulation;
class PositionalTruss : public WeakForm{
protected:
    // Material Young modulus
    double fYoungModulus;
    // Material Area
    double fArea;
public:
    //Positional elasticity 2D constructor
    PositionalTruss(int matid, int dim, double young, double area);

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

    /// @brief Computes the material constitutive matrix (used in plasticity models)
    /// @return constitutive matrix
    MatrixDouble ConstitutiveMatrix();

    /// @brief Returns the material young modulus
    /// @return Young modulus
    double YoungModulus(){
        return fYoungModulus;
    }

};


#endif