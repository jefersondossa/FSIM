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

    void ComputeStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, VecDouble &Rhs) override;
    
    void ComputeError(VecDouble &errors) override{};
    
    void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs) override;
};


#endif
