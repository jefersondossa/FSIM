//------------------------------------------------------------------------------
//-----------------------------------ELEMENT------------------------------------
//------------------------------------------------------------------------------

#ifndef ELEMENTMIXED_H
#define ELEMENTMIXED_H

#include "Element.h"
#include "MixedCompMesh.h"

/// Defines the fluid element object and all the element information
class ElementMixed : public Element{
protected:
    std::vector<Element *> fSubElements;     

public:
    ElementMixed() : Element(){
    };

    ElementMixed(int64_t index, std::vector<Element *> elvector, MixedCompMesh* mesh, WeakForm *wf);

    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) override;
    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override;

};


#endif

