#ifndef TRANSIENTWF_H
#define TRANSIENTWF_H

#include "WeakForm.h"
#include "CompMesh.h"

class TransientWeakForm
{
private:

public:
    TransientWeakForm() = default;
    ~TransientWeakForm() = default;

    virtual void UpdateTimeDerivatives(CompMesh *cmesh) = 0;
};


#endif