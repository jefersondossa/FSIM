#ifndef L2Projection_H
#define L2Projection_H

#include "WeakForm.h"

/// @brief Implements the class to enforce boundary conditions
class L2Projection : public WeakForm {
protected:
    // Boundary condition ID
    // 0 = Dirichlet (all directions)
    // 1 = Neumann (all directions)
    int BCType = 0;
    
    // Boundary condition value - Set into the stiffness matrix
    MatrixDouble BCVal1;
    
    // Boundary condition value - Set into the residual vector
    VecDouble    BCVal2;

public:
    /// @brief Boundary condition constructor
    /// @param matid physical tag
    /// @param dim problem dimension
    /// @param bctype BC type
    /// @param val1 value set into the stiffness matrix
    /// @param val2 value set into the residual vector
    L2Projection(int matid, int dim, int bctype, MatrixDouble &val1, VecDouble &val2);

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
    
    /// @brief Returns the BC value
    /// @return BC value
    VecDouble &BCValue(){return BCVal2;}
    
    /// @brief Post processing variable index - not used in this class
    /// @param name variable name
    /// @return variable index
    int VariableIndex(const std::string &name) const override{
        return 0;
    };

    /// @brief Returns the number of solution variables of a given post processing variable
    /// @param var solution variable's index
    /// @return number of solution variables
    int NSolutionVariables(int var) const override{
        return 0;
    };
};


#endif