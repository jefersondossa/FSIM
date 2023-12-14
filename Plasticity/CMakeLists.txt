
file(GLOB plasticity 
    PlasticityModel.h 
    PlasticityModel.cpp 
    LinearHardening.h
    LinearHardening.cpp
    MohrCoulomb.h 
    MohrCoulomb.cpp
    VonMises.h
    VonMises.cpp
    Tresca.h
    Tresca.cpp
    PlasticStep.h
    PlasticStep.cpp
)

set(PLASTICITY ${plasticity} PARENT_SCOPE)


