
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------------------ELEMENT------------------------------------
//------------------------------------------------------------------------------

#ifndef ELEMENTTRANSIENT_H
#define ELEMENTTRANSIENT_H

#include "ElementT.h"

/// Defines the fluid element object and all the element information
template<class geoshape, class compshape>
class ElementTransient : public ElementT<geoshape,compshape>{
    
public:
    ElementTransient() : ElementT<geoshape,compshape>(){
    };

    ElementTransient(int64_t index, VecInt &connect, CompMesh* mesh, WeakForm *wf);

    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) override;
    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override;

};


#endif

