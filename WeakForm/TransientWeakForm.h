#ifndef TRANSIENTWF_H
#define TRANSIENTWF_H

#include "WeakForm.h"

class TransientWeakForm
{
private:
    /* data */
public:
    TransientWeakForm() = default;
    ~TransientWeakForm() = default;

    virtual void UpdateTimeDerivatives() = 0;
};


#endif