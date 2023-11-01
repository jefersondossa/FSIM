#ifndef ELEMENTELASTICITYPOSITIONAL2D_H
#define ELEMENTELASTICITYPOSITIONAL2D_H

#include "ElementT.h"

template <class tshape>
class ElElasticityPositional2D : public ElementT<tshape>{
public:
    ElElasticityPositional2D(int index, VecInt &connect, CompMesh* mesh) : ElementT<tshape>(index,connect,mesh){};

    void ComputeStiffness(int &index, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, VecDouble &Rhs) override;
    
    void ComputeError(VecDouble &errors) override;
    
    void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs) override;

};


#endif