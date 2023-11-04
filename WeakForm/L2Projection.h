#ifndef L2Projection_H
#define L2Projection_H

#include "WeakForm.h"

class L2Projection : public WeakForm{
private:
    int fDimension;

    // Boundary condition ID
    int BCType = 0;
    
    // First value of boundary condition
    MatrixDouble BCVal1;
    
    // Second value of boundary condition
    VecDouble    BCVal2;

public:
    L2Projection(int matid, int dim, int bctype, MatrixDouble &val1, VecDouble &val2) : WeakForm(), fDimension(dim) {
        this->fMatId = matid;
        fNState = val2.size();
        BCType = bctype;
        BCVal1 = val1;
        BCVal2 = val2;
    };

    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override;
    
};


#endif