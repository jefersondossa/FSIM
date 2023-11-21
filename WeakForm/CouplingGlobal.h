#ifndef CouplingGlobal_H
#define CouplingGlobal_H

#include "WeakForm.h"
#include "CompMesh.h"

/// @brief Implements the Arlequin global coupling weak form
class CouplingGlobal : public WeakForm
{
private:
    // Global element index
    int64_t fGlobalIndex;
    // Global element computational mesh
    CompMesh * fGlobalMesh;
    // Global element adim coord
    MatrixDouble fGlobalXsi;
    // Global element correspondence to local mesh
    VecDouble fGlobalElemCorresp;
    // K0 and K1 in the classic Arlequin formulation (corresponding to L2 and H1(energy) portions)
    double fK0;
    double fK1;

public:
    // Global coupling Arlequin constructor
    CouplingGlobal(int dim, int64_t globindex, CompMesh* meshlocal, double k0, double k1);

    /// @brief Returns the element global index
    /// @return element global index
    int64_t &GetGlobalIndex(){
        return fGlobalIndex;
    }
    /// @brief Sets the global element parametric coordinates with respect to the local element
    /// @param globxsi parametric coordinates
    void SetGlobalXsi(MatrixDouble &globxsi){
        fGlobalXsi = globxsi;
    }

    /// @brief Sets the global element correspondence with respect to the local mesh. The vector stores the values of each integration point
    /// @param elemcorr local element correspondence
    void SetGlobalElemCorresp(VecDouble &elemcorr){
        fGlobalElemCorresp = elemcorr;
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
    void ArlequinStabStiffness(int &index, MatrixDouble &dphi_dx, VecDouble &phiGlobal, MatrixDouble &dphi_dxGlobal, double &weight_, double &djac_, std::vector<MatrixDouble> &Stiffness);

    /// @brief Computes the stabilized arlequin resudial vector 
    void ArlequinStabResidual(int &index, MatrixDouble &dphi_dx, VecDouble &phiGlobal, MatrixDouble &dphi_dxGlobal, double &weight_, double &djac_, std::vector<VecDouble> &Rhs);
};


#endif
