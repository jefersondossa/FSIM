#ifndef ElCouplingLocal_H
#define ElCouplingLocal_H

#include "ElementT.h"

template <class tshape>
class ElCouplingLocal : public ElementT<tshape>
{
private:
    int64_t fLocalIndex;
    std::vector<CompMesh *> fMeshVector;
    double tARLQ_ = 0.;

public:
    ElCouplingLocal(int index, int64_t fineindex, std::vector<CompMesh*> &meshvec):ElementT<tshape>(){
        this->Index()= index;
        this->fMesh = meshvec[2];
        fLocalIndex = fineindex;
        fMeshVector = meshvec;
        Element::nLocDOF = tshape::NElNodes * tshape::Dimension;
    };

    int64_t &GetLocalIndex(){return fLocalIndex;}

    void ComputeStiffness(int &index, std::vector<MatrixDouble> &Stiffness) override;
    
    void ComputeResidual(int &index, std::vector<VecDouble> &Rhs) override;
    
    void ComputeError(VecDouble &errors) override{};
    
    void ApplyBC(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override;
    
    void ArlequinStabStiffness(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, std::vector<MatrixDouble> &Stiffness);
    void ArlequinStabResidual(int &index, MatrixDouble &dphi_dx, double &weight_, double &djac_, std::vector<VecDouble> &Rhs);
};


#endif
