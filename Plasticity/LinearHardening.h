#ifndef LINEARHARDENING_H
#define LINEARHARDENING_H

#include "PlasticityModel.h"

class LinearHardening : public PlasticityModel
{
public:
    LinearHardening(int dim, double hardModulus);

};



#endif