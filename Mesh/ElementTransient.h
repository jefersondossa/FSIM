
//------------------------------------------------------------------------------
//-----------------------------------ELEMENT------------------------------------
//------------------------------------------------------------------------------

#ifndef ELEMENTTRANSIENT_H
#define ELEMENTTRANSIENT_H

#include "ElementT.h"

/// Defines the fluid element object and all the element information
template<class compshape>
class ElementTransient : public ElementT<compshape>{
    
public:
    ElementTransient() : ElementT<compshape>(){
    };

    ElementTransient(int index, GeoElement* gel, CompMesh* mesh, WeakForm *wf);

    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) override;
    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override;

};


#endif

