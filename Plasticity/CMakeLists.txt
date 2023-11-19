
file(GLOB plasticity 
    PlasticityModel.h 
    PlasticityModel.cpp 
    LinearHardening.h
    LinearHardening.cpp
)

set(PLASTICITY ${plasticity} PARENT_SCOPE)


