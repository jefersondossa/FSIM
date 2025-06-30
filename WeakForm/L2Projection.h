#ifndef L2Projection_H
#define L2Projection_H

#include "WeakForm.h"
#include "TransientWeakForm.h"

enum class BoundaryConditionType
{
    kDirichlet = 0,
    kNeumann = 1,
    kDirectionalHomogeneousDirichlet = 3,
    kDirectionalNonHomogeneousDirichlet = 4,
};

/// @brief Implements the class to enforce boundary conditions
class L2Projection : public WeakForm, public TransientWeakForm{
protected:
    // Boundary condition ID
    // 0 = Dirichlet (all directions)
    // 1 = Neumann (all directions)
    // 3 = Directional Homogeneous Dirichlet - apply Dirichlet BC in the non zero Val2 entry  
    // 4 = Directional Non-Homogeneous Dirichlet - apply Dirichlet BC in the non zero Val2 entry  
    BoundaryConditionType BCType{BoundaryConditionType::kDirichlet};
    
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
    L2Projection(int matid, int dim, BoundaryConditionType bctype, MatrixDouble &val1, VecDouble &val2);

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

    BoundaryConditionType GetBCType() const { return BCType; }
    
};


#endif