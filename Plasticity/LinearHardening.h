#ifndef LINEARHARDENING_H
#define LINEARHARDENING_H

#include "PlasticityModel.h"


class LinearHardening : public PlasticityModel
{
private:
    double fHardening;
    double fYield;
    MatrixDouble fConstitutiveMatrix;
    double fYoungModulus;

public:
    LinearHardening(WeakForm *elast, double hardModulus, double yield);

    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override;

    int VariableIndex(const std::string &name) const override;

    int NSolutionVariables(int var) const override;
    
    void Solution(IntPointData &data, int var, VecDouble &Sol) override;

    void ComputePlasticStrain(IntPointData &data, MatrixDouble &plasticstrain, MatrixDouble &totalstrain) override;
};



#endif