#ifndef ELEMENTNAVIERSTOKES_H
#define ELEMENTNAVIERSTOKES_H

#include "Stokes.h"

class NavierStokes : public Stokes{
public:
    NavierStokes(int matid, int dim, double density, double viscosity) : Stokes(matid,dim,density,viscosity) {
    
    };

    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override;

};


#endif