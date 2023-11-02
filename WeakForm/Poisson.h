#ifndef POISSON_H
#define POISSON_H

#include "WeakForm.h"

class Poisson : public WeakForm{
private:
    double  fScale = 1.;
    int     fDimension = 0;

public:
    Poisson(int matid, int dim, int nState = 1) : WeakForm(), fDimension(dim) {
        this->fMatId = matid;
        fNState = nState;
    };

    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override;
    
    void ApplyBC(std::vector<Node*> nodevec, VecInt &connect, MatrixDouble &Stiffness, VecDouble &Rhs) override;

    /// Sets the mesh moving weighting parameter for solving the Laplace problem
    /// @param double parameter value
    void SetScale(double &value) {fScale = value;};

    /// Gets the mesh moving weighting parameter
    /// @return mesh moving weighting parameter
    double &GetScale(){return fScale;};

    

    
};


#endif