#ifndef GlobalLocalEnrichment_H
#define GlobalLocalEnrichment_H

#include "WeakForm.h"


class GlobalLocalEnrichment : public WeakForm {
protected:
    // Young modulus
    REAL       fYoungModulus;
    // Poisson's ratio
    REAL       fPoissonRatio;
    // Thickness
    REAL       fThickness;
    // Problem consitutive matrix
    MatrixDouble fConstitutiveMatrix;
    // true if plane stress, false if plane strain
    bool         fPlaneStress;

public:
    GlobalLocalEnrichment(int matid, int dimension, REAL young, REAL poisson, bool planes = true, REAL thick=1);

    /// @brief Overloads the weak form stiffness matrix computation in the case more than one contribution is provided
    /// @param index integration point index
    /// @param data integration point data
    /// @param Stiffness vector of stiffness matrices
    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override{
        PanicButton();
    };
    void ComputeStiffness(int &index, IntPointData &localdata, IntPointData &globaldata, MatrixDouble &Stiffness);
    
    /// @brief Returns the weak form residual vector. It should never be called in this class
    /// @param index integration point index
    /// @param data integration point data
    /// @param Rhs residual vector
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override{
        PanicButton();
    };
    void ComputeResidual(int &index, IntPointData &localdata, IntPointData &globaldata, VecDouble &Rhs);
    
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