
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

set(WEAKFORM ${weakform} ${STOKES} ${NAVIERSTOKES} ${POISSON} ${ELASTICITY} ${ARLEQUIN} PARENT_SCOPE)


