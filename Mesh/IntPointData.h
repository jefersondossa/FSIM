#ifndef INTPOINTDATA
#define INTPOINTDATA

#include "DataTypes.h"


class IntPointData {

public:
    /// @brief Jacobian transform related to the initial mesh
    MatrixDouble fA0;
    /// @brief Inverse Jacobian transform related to the initial mesh
    MatrixDouble fA0Inv;
    /// @brief Element local axes related to the initial mesh
    MatrixDouble fAxes0;
   
    /// @brief Jacobian transform related to the current mesh
    MatrixDouble fA1;
    /// @brief Element local axes related to the current mesh
    MatrixDouble fAxes1;

    /// @brief Initial element Jacobian determinant
    double fJacA0;
    /// @brief Current element Jacobian determinant
    double fJacA1;

    /// @brief Integration weight
    double fWeight;
    /// @brief Integration point coordinates in local axes
    VecDouble fAdimCoord;

    /// @brief Integration point in the initial mesh
    VecDouble fX;
    /// @brief Integration point in the current mesh
    VecDouble fX1;
    
    /// @brief Shape functions
    VecDouble    fPhi;
    /// @brief Shape functions derivatives related to the local axes
    MatrixDouble fDPhi;
    /// @brief Shape function derivatives related to the initial mesh
    MatrixDouble fDPhiX0;
    /// @brief Shape function derivatives related to the current mesh
    MatrixDouble fDPhiX1;

    /// @brief Arlequin Weight function
    VecDouble fWeightFunction;
    /// @brief Previous time step Arlequin weight function
    VecDouble fPrevWeightFunction;
    /// @brief Signaled distance
    VecDouble fDistFunction;

    /// @brief Solution
    VecDouble fSol;
    /// @brief Solution gradient
    MatrixDouble fDSolDx;
    /// @brief Solution first time derivative
    VecDouble fDSolDt;
    /// @brief Solution second time derivative
    VecDouble fDSolDDt;
    /// @brief Previous time/load step Solution
    VecDouble fSolPrev;

    /// @brief Flags implemented to compute or not the solution and its derivatives
    bool fNeedsSol = false;
    bool fNeedsDSol = false;
    bool fNeedsTimeDerivatives = false;

    VecDouble fPlasticMultiplier;
    VecDouble fPlasticStrain;
    VecDouble fEffectiveStress;

};


#endif