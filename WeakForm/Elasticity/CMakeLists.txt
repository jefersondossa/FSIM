
file(GLOB elasticity 
    Elasticity2D.h
    Elasticity2D.cpp
    Elasticity3D.h
    Elasticity3D.cpp
    ElasticityPositional2D.h
    ElasticityPositional2D.cpp
    ElasticTruss.h 
    ElasticTruss.cpp
    PositionalTruss.h
    PositionalTruss.cpp
    TransientElasticity2D.h
    TransientElasticity2D.cpp
    TransientElasticTruss.h
    TransientElasticTruss.cpp
    TransientPositionalTruss.h
    TransientPositionalTruss.cpp
    TransientPositional2D.h
    TransientPositional2D.cpp
    LinearBeam.h
    LinearBeam.cpp
    LinearFrame.h
    LinearFrame.cpp
    PositionalFrame2D.h
    PositionalFrame2D.cpp
)

set(ELASTICITY ${elasticity} PARENT_SCOPE)


