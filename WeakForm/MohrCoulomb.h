#ifndef MOHRCOULOMB_H
#define MOHRCOULOMB_H

#include "Elasticity2D.h"

class MohrCoulomb : public Elasticity2D{
protected:
    double fInternalFriction;
    double fCohesion;

public:
    MohrCoulomb(int matid, double young, double poisson, double friction, double cohesion, bool planes = true) : 
                Elasticity2D(matid,young,poisson,planes), fInternalFriction(friction), fCohesion(cohesion) {
    };

    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override;
};


#endif