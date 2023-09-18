#ifndef ELEMENTPOISSON_H
#define ELEMENTPOISSON_H

#include "Element.h"

class ElPoisson : public Element{

public:
    ElPoisson(int index, VecInt &connect, CompMesh* mesh) : Element(index,connect,mesh){};

    void ComputeStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, VecDouble &Rhs) override;
    
    void ComputeError(VecDouble &errors) override;
    
    void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs) override;


};


#endif