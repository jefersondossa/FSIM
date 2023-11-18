#ifndef ELASTICTRUSS_H
#define ELASTICTRUSS_H

#include "WeakForm.h"

class ElasticTruss : public WeakForm{
protected:
    // enum 
    bool         fPlaneStress;
    double       fYoungModulus;
    double       fArea;

public:    
    ElasticTruss(int matid, int dim, double young, double area);

    void ComputeStiffness(int &index, IntPointData &data, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, VecDouble &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override;

    int VariableIndex(const std::string &name) const override;

    int NSolutionVariables(int var) const override;
    
    void Solution(IntPointData &data, int var, VecDouble &Sol) override;
};


#endif