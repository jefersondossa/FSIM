#ifndef ESO_HPP
#define ESO_HPP
#include "DenseEigen.h"

class CompMesh;
class LinearAnalysis;

void RunEso(CompMesh& model, LinearAnalysis& an, REAL target_final_vol = 0.7);

#endif