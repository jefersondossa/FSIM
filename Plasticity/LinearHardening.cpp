#include "LinearHardening.h"

LinearHardening::LinearHardening(int dim, double hardModulus) : PlasticityModel(){
    fHardening = hardModulus;
    fTotalStrain.resize(dim,dim);
    fTotalStrain.setZero();
    fPlastStrain.resize(dim,dim);
    fPlastStrain.setZero();
}