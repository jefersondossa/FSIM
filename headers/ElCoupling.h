#ifndef ELCOUPLING_H
#define ELCOUPLING_H

#include "Element.h"

class ElCoupling : public Element
{
private:
    int64_t fLocalIndex;
    int64_t fGlobalIndex;
    std::vector<CompMesh *> fMeshVector;

public:
    ElCoupling(int index, int64_t fineindex, std::vector<CompMesh*> &meshvec):Element(){
        this->Index()= index;
        SetMesh(meshvec[2]);
        fLocalIndex = fineindex;
        fMeshVector = meshvec;
    };

    int64_t &GetLocalIndex(){return fLocalIndex;}
    int64_t &GetGlobalIndex(){return fGlobalIndex;}
    void SetGlobalIndex(int64_t gindex){fGlobalIndex = gindex;}

    void ComputeStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, MatrixDouble &Stiffness) override;
    
    void ComputeResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, VecDouble &Rhs) override{};
    
    void ComputeError(VecDouble &errors) override{};
    
    void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs) override{};
};


#endif
