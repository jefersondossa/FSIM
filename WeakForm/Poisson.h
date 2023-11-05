#ifndef POISSON_H
#define POISSON_H

#include "WeakForm.h"

class Poisson : public WeakForm{
private:
    double  fScale = 1.;
    enum SolutionVariables{ESolution, EDerivative};
public:
    Poisson(int matid, int dim, int nState = 1) : WeakForm() {
        this->fMatId = matid;
        this->fDimension = dim;
        fNState = nState;
    };

    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override;
    
    /// Sets the mesh moving weighting parameter for solving the Laplace problem
    /// @param double parameter value
    void SetScale(double &value) {fScale = value;};

    /// Gets the mesh moving weighting parameter
    /// @return mesh moving weighting parameter
    double &GetScale(){return fScale;};

    

    
};


#endif