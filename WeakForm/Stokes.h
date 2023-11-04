#ifndef ELEMENTSTOKES_H
#define ELEMENTSTOKES_H

#include "WeakForm.h"

class Stokes : public WeakForm{
protected:
    double fViscosity;
    double fDensity;
    void GetStabilizationParameter(int &index, IntPointData &data);
    double tPSPG_, tSUPG_, tLSIC_; 

public:
    Stokes(int matid, int dim, double density, double viscosity) : WeakForm() {
        this->fMatId = matid;
        fDimension = dim;
        fNState = dim+1;
        fViscosity = viscosity;
        fDensity = density;
    };

    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override;

};


#endif