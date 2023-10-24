#ifndef ELEMENTSTOKES_H
#define ELEMENTSTOKES_H

#include "Element.h"

class ElStokes : public Element{
public:
    ElStokes(int index, VecInt &connect, CompMesh* mesh) : Element(index,connect,mesh){};

    void ComputeStiffness(int &index, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, VecDouble &Rhs) override;
    
    void ComputeError(VecDouble &errors) override;
    
    void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs) override;

// private:
    void GetStabilizationParameter(int &index, double &tSUPG_, double &tPSPG_, double &tLSIC_, MatrixDouble &dphi_dx);
    double tPSPG_, tSUPG_, tLSIC_; 

};


#endif