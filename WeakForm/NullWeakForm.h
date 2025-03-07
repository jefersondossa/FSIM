#ifndef NullWF_H
#define NullWF_H

#include "WeakForm.h"


/// @brief Implements the class to enforce boundary conditions
class NullWeakForm : public WeakForm {
public:
    /// @brief Boundary condition constructor
    /// @param matid physical tag
    /// @param dim problem dimension
    /// @param bctype BC type
    /// @param val1 value set into the stiffness matrix
    /// @param val2 value set into the residual vector
    NullWeakForm(int matid, int nstate);

    /// @brief Overloads the weak form stiffness matrix computation in the case more than one contribution is provided
    /// @param index integration point index
    /// @param data integration point data
    /// @param Stiffness vector of stiffness matrices
    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override{
        // do nothing
        return;
    };
    
    /// @brief Returns the weak form residual vector. It should never be called in this class
    /// @param index integration point index
    /// @param data integration point data
    /// @param Rhs residual vector
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override{
        //do nothing
        return;
    };

};


#endif