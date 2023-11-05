#ifndef CouplingGlobal_H
#define CouplingGlobal_H

#include "WeakForm.h"
#include "CompMesh.h"

class CouplingGlobal : public WeakForm
{
private:
    int64_t fGlobalIndex;
    CompMesh * fGlobalMesh;
    MatrixDouble fGlobalXsi;
    VecDouble fGlobalElemCorresp;
    double fK0;
    double fK1;

public:
    CouplingGlobal(int dim, int64_t globindex, CompMesh* meshlocal, double k0, double k1):WeakForm(){
        fGlobalMesh = meshlocal;
        fDimension = dim;
        fK0 = k0;
        fK1 = k1;
        fGlobalIndex = globindex;
        
    };

    int64_t &GetGlobalIndex(){
        return fGlobalIndex;
    }
    void SetGlobalXsi(MatrixDouble &globxsi){
        fGlobalXsi = globxsi;
    }
    void SetGlobalElemCorresp(VecDouble &elemcorr){
        fGlobalElemCorresp = elemcorr;
    }

    void ComputeStiffness(int &index, IntPointData &data, std::vector<MatrixDouble> &Stiffness) override;
    
    void ComputeResidual(int &index, IntPointData &data, std::vector<VecDouble> &Rhs) override;
    
    void ComputeError(IntPointData &data, VecDouble &errors) override{};
    
    void ArlequinStabStiffness(int &index, MatrixDouble &dphi_dx, VecDouble &phiGlobal, MatrixDouble &dphi_dxGlobal, double &weight_, double &djac_, std::vector<MatrixDouble> &Stiffness);
    void ArlequinStabResidual(int &index, MatrixDouble &dphi_dx, VecDouble &phiGlobal, MatrixDouble &dphi_dxGlobal, double &weight_, double &djac_, std::vector<VecDouble> &Rhs);
};


#endif
