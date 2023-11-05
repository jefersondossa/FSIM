
file(GLOB weakform 
    CouplingGlobal.h
    CouplingGlobal.cpp
    CouplingLocal.h
    CouplingLocal.cpp
    Elasticity2D.h
    Elasticity2D.cpp
    ElasticityPositional2D.h
    ElasticityPositional2D.cpp
    NavierStokes.h
    NavierStokes.cpp
    Stokes.h
    Stokes.cpp
    WeakForm.h
    WeakForm.cpp
    Poisson.h
    Poisson.cpp
    L2Projection.h 
    L2Projection.cpp
)

set(WEAKFORM ${weakform} PARENT_SCOPE)


