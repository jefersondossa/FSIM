#ifndef LagrangeMultiplier_H
#define LagrangeMultiplier_H

#include "WeakForm.h"


/// @brief Implements the class to enforce boundary conditions
class LagrangeMultiplier : public WeakForm {
protected:
    double fScale = 1.0;

public:
    /// @brief Boundary condition constructor
    /// @param matid physical tag
    /// @param dim problem dimension
    /// @param bctype BC type
    /// @param val1 value set into the stiffness matrix
    /// @param val2 value set into the residual vector
    LagrangeMultiplier(int matid, int nstate, double scale = 1.0);

    /// @brief Overloads the weak form stiffness matrix computation in the case more than one contribution is provided
    /// @param index integration point index
    /// @param data integration point data
    /// @param Stiffness vector of stiffness matrices
    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override{
        PanicButton();
    };
    void ComputeStiffness(int &index, IntPointData &leftdata, IntPointData &rightdata, MatrixDouble &Stiffness);
    
    /// @brief Returns the weak form residual vector. It should never be called in this class
    /// @param index integration point index
    /// @param data integration point data
    /// @param Rhs residual vector
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override{
        PanicButton();
    };
    void ComputeResidual(int &index, IntPointData &leftdata, IntPointData &rightdata, VecDouble &Rhs);
    
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

};


#endif