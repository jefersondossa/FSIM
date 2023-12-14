#include "PlasticStep.h"

template<class TPlastic>
PlasticStep<TPlastic>::PlasticStep(/* args */)
{
}

template<class TPlastic>
PlasticStep<TPlastic>::~PlasticStep()
{
}

#include "MohrCoulomb.h"
#include "Tresca.h"
#include "VonMises.h"

template class PlasticStep<MohrCoulomb>;
template class PlasticStep<Tresca>;
template class PlasticStep<VonMises>;