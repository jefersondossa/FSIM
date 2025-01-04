#ifndef ESO_HPP
#define ESO_HPP

class CompMesh;
class LinearAnalysis;

void RunEso(CompMesh& model, LinearAnalysis& an, double target_final_vol = 0.7);

#endif