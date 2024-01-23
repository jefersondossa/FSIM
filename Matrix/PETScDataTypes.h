#ifndef PETSC_DATA_TYPES
#define PETSC_DATA_TYPES

#include <cmath>
#include <stdio.h>
#include <vector>
#include <functional>

#ifdef HAS_PETSC

#include <petscksp.h>

/// @brief Petsc Matrix
typedef Mat PETScMat;
/// @brief Petsc vector
typedef Vec PETScVec;

#endif

#endif