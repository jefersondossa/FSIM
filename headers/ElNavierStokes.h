#ifndef ELEMENTNAVIERSTOKES_H
#define ELEMENTNAVIERSTOKES_H

#include "ElStokes.h"

template <class tshape>
class ElNavierStokes : public ElStokes<tshape>{
public:
    ElNavierStokes(int index, VecInt &connect, CompMesh* mesh) : ElStokes<tshape>(index,connect,mesh){};

    void ComputeStiffness(int &index, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, VecDouble &Rhs) override;
    
    void ComputeError(VecDouble &errors) override;
};


#endif