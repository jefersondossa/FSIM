#ifndef PLASTICITYMODEL_H
#define PLASTICITYMODEL_H

#include "DataTypes.h"

class PlasticityModel
{
protected: 
    MatrixDouble fTotalStrain;
    MatrixDouble fPlastStrain;

    double fHardening;

    bool fIsPlastic = false;

public:
    PlasticityModel() = default;
    ~PlasticityModel(){};
};



#endif