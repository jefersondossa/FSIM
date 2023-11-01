#ifndef ELEMENTSTOKES_H
#define ELEMENTSTOKES_H

#include "ElementT.h"

template <class tshape>
class ElStokes : public ElementT<tshape>{
public:
    ElStokes(int index, VecInt &connect, CompMesh* mesh) : ElementT<tshape>(index,connect,mesh){};

    void ComputeStiffness(int &index, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, VecDouble &Rhs) override;
    
    void ComputeError(VecDouble &errors) override;
    
    void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs) override;

protected:
    void GetStabilizationParameter(int &index, double &tSUPG_, double &tPSPG_, double &tLSIC_, MatrixDouble &dphi_dx);
    
    double tPSPG_, tSUPG_, tLSIC_; 
};


#endif