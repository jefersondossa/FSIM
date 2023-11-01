#ifndef ElCouplingGlobal_H
#define ElCouplingGlobal_H

#include "ElementT.h"

template <class tshape>
class ElCouplingGlobal : public ElementT<tshape>
{
private:
    int64_t fGlobalIndex;
    std::vector<CompMesh *> fMeshVector;
    MatrixDouble fGlobalXsi;
    VecDouble fGlobalElemCorresp;

public:
    ElCouplingGlobal(int index, int64_t globindex, std::vector<CompMesh*> &meshvec):ElementT<tshape>(){
        this->Index()= index;
        this->fMesh = meshvec[2];
        fGlobalIndex = globindex;
        fMeshVector = meshvec;
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

    void ComputeStiffness(int &index, std::vector<MatrixDouble> &Stiffness) override;
    
    void ComputeResidual(int &index, std::vector<VecDouble> &Rhs) override;
    
    void ComputeError(VecDouble &errors) override{};
    
    void ApplyBC(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override;

    void ArlequinStabStiffness(int &index, MatrixDouble &dphi_dx, VecDouble &phiGlobal, MatrixDouble &dphi_dxGlobal, double &weight_, double &djac_, std::vector<MatrixDouble> &Stiffness);
    void ArlequinStabResidual(int &index, MatrixDouble &dphi_dx, VecDouble &phiGlobal, MatrixDouble &dphi_dxGlobal, double &weight_, double &djac_, std::vector<VecDouble> &Rhs);
};


#endif
