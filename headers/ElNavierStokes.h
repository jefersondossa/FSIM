#ifndef ELEMENTNAVIERSTOKES_H
#define ELEMENTNAVIERSTOKES_H

#include "ElStokes.h"

class ElNavierStokes : public ElStokes{
public:
    ElNavierStokes(int index, VecInt &connect, CompMesh* mesh) : ElStokes(index,connect,mesh){};

    void ComputeStiffness(int &index, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, VecDouble &Rhs) override;
    
    void ComputeError(VecDouble &errors) override;
};


#endif