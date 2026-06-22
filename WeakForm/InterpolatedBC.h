#ifndef InterpolatedBC_H
#define InterpolatedBC_H

#include "WeakForm.h"
#include "TransientWeakForm.h"
#include "L2Projection.h"

/// @brief Implements the class to enforce boundary conditions
class InterpolatedBC : public WeakForm, public TransientWeakForm{
protected:
    // Boundary condition ID
    // 0 = Dirichlet (all directions)
    // 1 = Neumann (all directions)
    // 3 = Directional Homogeneous Dirichlet - apply Dirichlet BC in the non zero Val2 entry  
    // 4 = Directional Non-Homogeneous Dirichlet - apply Dirichlet BC in the non zero Val2 entry  
    BoundaryConditionType BCType{BoundaryConditionType::kDirichlet};
    
    std::map<int,int> *fGlobalElementCorrespondence;
    std::map<int, MatrixDouble> *fGlobalNodeCorrespondence;
    CompMesh *fGlobalMesh;

public:
    /// @brief Boundary condition constructor
    /// @param matid physical tag
    /// @param dim problem dimension
    /// @param bctype BC type
    InterpolatedBC(int matid, int dim, int nstate, BoundaryConditionType bctype, 
        std::map<int,int> *globalElementCorrespondence,
        std::map<int, MatrixDouble> *globalNodeCorrespondence,
        CompMesh *fGlobalMesh);

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

    BoundaryConditionType GetBCType() const { return BCType; }
    
};


#endif