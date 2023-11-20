#ifndef MOHRCOULOMB_H
#define MOHRCOULOMB_H

#include "PlasticityModel.h"


class MohrCoulomb : public PlasticityModel
{
private:
    double fHardening;
    double fYield;
    MatrixDouble fConstitutiveMatrix;
    double fInternalFriction;
    double fCohesion;
    double fYoungModulus;
    MatrixDouble fPlasticStrain;
    MatrixDouble fTotalStrain;

public:
    LinearHardening(WeakForm *elast, double hardModulus, double intfriction);

    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override;

    int VariableIndex(const std::string &name) const override;

    int NSolutionVariables(int var) const override;
    
    void Solution(IntPointData &data, int var, VecDouble &Sol) override;

    void ComputePlasticStrain(IntPointData &data, MatrixDouble &plasticstrain, MatrixDouble &totalstrain) override;
};



#endif