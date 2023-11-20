
file(GLOB plasticity 
    PlasticityModel.h 
    PlasticityModel.cpp 
    LinearHardening.h
    LinearHardening.cpp
    MohrCoulomb.h 
    MohrCoulomb.cpp
)

set(PLASTICITY ${plasticity} PARENT_SCOPE)


