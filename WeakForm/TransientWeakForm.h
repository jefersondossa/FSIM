#ifndef TRANSIENTWF_H
#define TRANSIENTWF_H

#include "WeakForm.h"
#include "CompMesh.h"

class TransientWeakForm
{
public:
    enum TimeIntegScheme {ENewmark, EGeneralizedAlpha};
protected:

    TimeIntegScheme fIntegScheme;

    // Time Step
    REAL fTimeStep = 0;
    int fTimeInstant = 0;

    // Lambda function returning the exact solution and its gradient
    std::function<void (const VecDouble &coord, REAL &time, VecDouble &u, MatrixDouble &gradU)> fExactSolTransient = 0; 
    // Lambda function returning the forcing function
    std::function<void (const VecDouble &coord, REAL &time, VecDouble &force)> fForceFunctionTransient = 0; 

public:
    TransientWeakForm() = default;
    ~TransientWeakForm() = default;

    virtual void UpdateTimeDerivatives(CompMesh *cmesh){
        //Should never be called in this class
        PanicButton();
    };

    void SetTimeStep(REAL dTime){
        fTimeStep = dTime;
    }
    void SetTimeInstant(int instant){
        fTimeInstant = instant;
    }

    /// @brief Sets the exact solution
    /// @param exSol Lambda function to compute the exact solution and its gradient for a given coordinate
    void SetExactSolutionTransient(std::function<void (const VecDouble &coord, REAL &time, VecDouble &u, MatrixDouble &gradU)> exSol){
        fExactSolTransient = exSol;
    }

    /// @brief Gets the exact solution Lambda function
    /// @return Exact solution Lambda function
    std::function<void (const VecDouble &coord, REAL &time, VecDouble &u, MatrixDouble &gradU)> &GetExactSolutionTransient(){
        return fExactSolTransient;
    }

    /// @brief Sets the source term Lambda function
    /// @param ffunction Lambda function to compute the source term for given coordinate.
    void SetForcingFunctionTransient(std::function<void (const VecDouble &coord, REAL &time, VecDouble &force)> ffunction){
        fForceFunctionTransient = ffunction;
    }

    /// @brief Gets the source term Lambda function
    /// @return Source term Lambda function
    std::function<void (const VecDouble &coord, REAL &time, VecDouble &force)> &GetForcingFunctionTransient(){
        return fForceFunctionTransient;
    }
};


#endif