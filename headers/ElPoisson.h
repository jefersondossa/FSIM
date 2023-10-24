#ifndef ELEMENTPOISSON_H
#define ELEMENTPOISSON_H

#include "Element.h"

class ElPoisson : public Element{
private:
    double  fMeshMovingParameter = 1.;

public:
    ElPoisson(int index, VecInt &connect, CompMesh* mesh) : Element(index,connect,mesh){};

    void ComputeStiffness(int &index, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, VecDouble &Rhs) override;
    
    void ComputeError(VecDouble &errors) override;
    
    void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs) override;

    /// Sets the mesh moving weighting parameter for solving the Laplace problem
    /// @param double parameter value
    void setMeshMovingParameter(double &value) {fMeshMovingParameter = value;};

    /// Gets the mesh moving weighting parameter
    /// @return mesh moving weighting parameter
    double &getMeshMovingParameter(){return fMeshMovingParameter;};

};


#endif