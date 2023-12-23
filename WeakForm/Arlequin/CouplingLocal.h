#ifndef ElCouplingLocal_H
#define ElCouplingLocal_H

#include "WeakForm.h"
#include "CompMesh.h"

/// @brief Implements the Arlequin local coupling operator
class CouplingLocal : public WeakForm {
protected:
    // Local element index
    int64_t fLocalIndex;
    // Local computational mesh
    CompMesh * fLocalMesh;
    // Arlequin stabilization parameter
    double tARLQ_ = 0.;
    // K0 and K1 in the classic Arlequin formulation (corresponding to L2 and H1(energy) portions)
    double fK0;
    double fK1;

public:
    // Local coupling constructor
    CouplingLocal(int dim, int64_t fineindex, CompMesh* meshlocal, double k0, double k1);

    /// @brief Returns the element local index
    /// @return local mesh element index
    int64_t &GetLocalIndex(){
        return fLocalIndex;
    }

    /// @brief Overloads the weak form stiffness matrix computation in the case more than one contribution is provided
    /// @param index integration point index
    /// @param data integration point data
    /// @param Stiffness vector of stiffness matrices
    void ComputeStiffness(int &index, IntPointData &data, std::vector<MatrixDouble> &Stiffness) override;
    
    /// @brief Returns the weak form residual vector. It should never be called in this class
    /// @param index integration point index
    /// @param data integration point data
    /// @param Rhs residual vector
    void ComputeResidual(int &index, IntPointData &data, std::vector<VecDouble> &Rhs) override;
    
    /// @brief Computes the element error. The exact solution shoul be provided.
    /// @param data integration point data
    /// @param errors vector storing all errors
    void ComputeError(IntPointData &data, VecDouble &errors) override {};
    
    /// @brief Computes the stabilized Arlequin stiffness matrix
    void ArlequinStabStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, std::vector<MatrixDouble> &Stiffness);
    
    /// @brief Computes the stabilized arlequin resudial vector 
    void ArlequinStabResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, std::vector<VecDouble> &Rhs);
};


#endif
