#ifndef ElCouplingLocal_H
#define ElCouplingLocal_H

#include "Element.h"

class ElCouplingLocal : public Element
{
private:
    int64_t fLocalIndex;
    std::vector<CompMesh *> fMeshVector;

public:
    ElCouplingLocal(int index, int64_t fineindex, std::vector<CompMesh*> &meshvec):Element(){
        this->Index()= index;
        SetMesh(meshvec[2]);
        fLocalIndex = fineindex;
        fMeshVector = meshvec;
    };

    int64_t &GetLocalIndex(){return fLocalIndex;}

    void ComputeStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, VecDouble &Rhs) override;
    
    void ComputeError(VecDouble &errors) override{};
    
    void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs) override;
};


#endif
