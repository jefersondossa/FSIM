#ifndef ELEMENTELASTICITY2D_H
#define ELEMENTELASTICITY2D_H

#include "ElementT.h"

template <class tshape>
class ElElasticity2D : public ElementT<tshape>{
public:
    ElElasticity2D(int index, VecInt &connect, CompMesh* mesh) : ElementT<tshape>(index,connect,mesh){
        Element::nLocDOF = tshape::NElNodes * tshape::Dimension;
    };

    void ComputeStiffness(int &index, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, VecDouble &Rhs) override;
    
    void ComputeError(VecDouble &errors) override;
    
    void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs) override;

};


#endif