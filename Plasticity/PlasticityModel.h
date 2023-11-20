#ifndef PLASTICITYMODEL_H
#define PLASTICITYMODEL_H

#include "DataTypes.h"
#include "WeakForm.h"
#include "L2Projection.h"

class PlasticityModel : public WeakForm
{
protected: 
    int fRealDimension;

    WeakForm *fElasticModel;

public:
    PlasticityModel() = default;
    PlasticityModel(WeakForm *elast);
    ~PlasticityModel(){};

    void BuildStressTensor(VecDouble &Sol, MatrixDouble &Stress);
    void ComputePrincipalStress(MatrixDouble &Stress, VecDouble &PrincipalS);

    virtual void ComputePlasticStrain(IntPointData &data, MatrixDouble &plasticstrain, MatrixDouble &totalstrain) = 0;

    int &RealDimension(){return fRealDimension;}
};



#endif