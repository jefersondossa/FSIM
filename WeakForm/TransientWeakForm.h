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
    double fTimeStep;

public:
    TransientWeakForm() = default;
    ~TransientWeakForm() = default;

    virtual void UpdateTimeDerivatives(CompMesh *cmesh) = 0;
};


#endif