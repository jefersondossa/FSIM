
//------------------------------------------------------------------------------
//-----------------------------------ELEMENT------------------------------------
//------------------------------------------------------------------------------

#ifndef ELEMENTWITHMEM_H
#define ELEMENTWITHMEM_H

#include "ElementT.h"
#include "PlasticityModel.h"

/// Defines the fluid element object and all the element information
template<class compshape>
class ElementWithMem : public ElementT<compshape>{
protected:
    MatrixDouble fElasticConstitutiveMatrix;

public:
    ElementWithMem() : ElementT<compshape>(){
    };

    ElementWithMem(int index, GeoElement* gel, CompMesh* mesh, WeakForm *wf);

    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) override;
    void ComputeElContribution(MatrixDouble &Stiffness) override;
    void ComputeElContribution(VecDouble &Rhs) override;
    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override;
    // void ComputeElContribution(std::vector<MatrixDouble> &Stiffness) override;
    // void ComputeElContribution(std::vector<VecDouble> &Rhs) override;

    void ComputeTrialStress(int &index,Tensor3D &ElasStress);

};


#endif

