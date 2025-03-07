
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

#ifndef ELEMENTLAGRANGEMULT_H
#define ELEMENTLAGRANGEMULT_H

#include "ElementT.h"

/// Defines the fluid element object and all the element information
template<class tshape>
class ElementLagrangeMultiplier : public ElementT<tshape>{
    Element* fLeftElement;
    Element* fRightElement;
    
public:
    ElementLagrangeMultiplier() : ElementT<tshape>(){
    };

    ElementLagrangeMultiplier(int64_t index, Element* leftEl, Element* rightEl, CompMesh* mesh, WeakForm *wf);

    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) override;

    void GetInterfaceData(Element* &leftel, Element* &rightel) override{
        leftel = fLeftElement;
        rightel = fRightElement;
    };

};


#endif

