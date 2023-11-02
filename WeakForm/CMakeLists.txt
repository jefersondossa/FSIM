
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
    ElPoisson.h
    ElPoisson.cpp
)

set(WEAKFORM ${weakform} PARENT_SCOPE)


