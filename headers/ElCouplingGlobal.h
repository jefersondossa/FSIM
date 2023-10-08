#ifndef ElCouplingGlobal_H
#define ElCouplingGlobal_H

#include "Element.h"

class ElCouplingGlobal : public Element
{
private:
    int64_t fGlobalIndex;
    std::vector<CompMesh *> fMeshVector;
    MatrixDouble fGlobalXsi;
    VecDouble fGlobalElemCorresp;

public:
    ElCouplingGlobal(int index, int64_t globindex, std::vector<CompMesh*> &meshvec):Element(){
        this->Index()= index;
        SetMesh(meshvec[2]);
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

    void ComputeStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, std::vector<MatrixDouble> &Stiffness) override;
    
    void ComputeResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, std::vector<VecDouble> &Rhs) override;
    
    void ComputeError(VecDouble &errors) override{};
    
    void ApplyBC(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override;

    void ArlequinStabStiffness(int &index, MatrixDouble &dphi_dx, VecDouble &phiGlobal, MatrixDouble &dphi_dxGlobal, double &weight_, double &djac_, std::vector<MatrixDouble> &Stiffness);
    void ArlequinStabResidual(int &index, MatrixDouble &dphi_dx, VecDouble &phiGlobal, MatrixDouble &dphi_dxGlobal, double &weight_, double &djac_, std::vector<VecDouble> &Rhs);
};


#endif
