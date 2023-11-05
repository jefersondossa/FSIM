#ifndef ElCouplingLocal_H
#define ElCouplingLocal_H

#include "WeakForm.h"
#include "CompMesh.h"

class CouplingLocal : public WeakForm
{
private:
    int64_t fLocalIndex;
    CompMesh * fLocalMesh;
    double tARLQ_ = 0.;
    double fK0;
    double fK1;

public:
    CouplingLocal(int dim, int64_t fineindex, CompMesh* meshlocal, double k0, double k1) : WeakForm(){
        fLocalMesh = meshlocal;
        fLocalIndex = fineindex;
        fDimension = dim;
        fK0 = k0;
        fK1 = k1;
    };

    int64_t &GetLocalIndex(){return fLocalIndex;}

    void ComputeStiffness(int &index, IntPointData &data, std::vector<MatrixDouble> &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, std::vector<VecDouble> &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override{};
        
    void ArlequinStabStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, std::vector<MatrixDouble> &Stiffness);
    void ArlequinStabResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, std::vector<VecDouble> &Rhs);
};


#endif
