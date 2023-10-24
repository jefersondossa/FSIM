#ifndef ELEMENTELASTICITYPOSITIONAL2D_H
#define ELEMENTELASTICITYPOSITIONAL2D_H

#include "Element.h"

class ElElasticityPositional2D : public Element{
public:
    ElElasticityPositional2D(int index, VecInt &connect, CompMesh* mesh) : Element(index,connect,mesh){};

    void ComputeStiffness(int &index, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, VecDouble &Rhs) override;
    
    void ComputeError(VecDouble &errors) override;
    
    void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs) override;

};


#endif