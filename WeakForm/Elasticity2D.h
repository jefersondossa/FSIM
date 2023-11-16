#ifndef ELASTICITY2D_H
#define ELASTICITY2D_H

#include "WeakForm.h"

class Elasticity2D : public WeakForm{
protected:
    bool fPlaneStress = true;
    double fYoungModulus;
    double fPoissonRatio;

public:
    Elasticity2D(int matid, double young, double poisson, bool planes = true) : WeakForm() {
        this->fMatId = matid;
        fDimension = 2;
        fNState = 2;
        fYoungModulus = young;
        fPoissonRatio = poisson;
        fPlaneStress = planes;
    };

    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override;
};


#endif