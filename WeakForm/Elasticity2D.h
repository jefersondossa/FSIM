#ifndef ELASTICITY2D_H
#define ELASTICITY2D_H

#include "WeakForm.h"

class Elasticity2D : public WeakForm{
protected:
    // enum 
    bool         fPlaneStress;
    double       fYoungModulus;
    double       fPoissonRatio;
    MatrixDouble fConstitutiveMatrix;

public:
    enum PostProcVar {EDisplacement, EStress, EStrain, EExactStress, EExactDisplacement, EExactStrain};
    
    Elasticity2D(int matid, double young, double poisson, bool planes = true);

    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override;

    int VariableIndex(const std::string &name) const override;

    int NSolutionVariables(int var) const override;
    
    void Solution(IntPointData &data, int var, VecDouble &Sol) override;
};


#endif