
file(GLOB weakform 
    ElCouplingGlobal.h
    ElCouplingGlobal.cpp
    ElCouplingLocal.h
    ElCouplingLocal.cpp
    ElElasticity2D.h
    ElElasticity2D.cpp
    ElElasticityPositional2D.h
    ElElasticityPositional2D.cpp
    ElNavierStokes.h
    ElNavierStokes.cpp
    ElStokes.h
    ElStokes.cpp
    WeakForm.h
    WeakForm.cpp
    Poisson.h
    Poisson.cpp
    L2Projection.h 
    L2Projection.cpp
)

set(WEAKFORM ${weakform} PARENT_SCOPE)


