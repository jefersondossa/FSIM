
file(GLOB weakform 
    WeakForm.h
    WeakForm.cpp
    L2Projection.h 
    L2Projection.cpp
    TransientWeakForm.h
    TransientWeakForm.cpp
)

add_subdirectory(Stokes)
add_subdirectory(NavierStokes)
add_subdirectory(Elasticity)
add_subdirectory(Arlequin)
add_subdirectory(Poisson)
add_subdirectory(PhaseField)

set(WEAKFORM ${weakform} ${STOKES} ${NAVIERSTOKES} ${POISSON} ${PHASE_FIELD} ${ELASTICITY} ${ARLEQUIN} PARENT_SCOPE)


