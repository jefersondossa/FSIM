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

    virtual Element *Clone() const;

    ElementMixed(int64_t index, std::vector<Element *> elvector, MixedCompMesh* mesh, WeakForm *wf);

    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) override;
    void ComputeElContribution(VecDouble &Rhs) override{
        PanicButton();
    };
    void ComputeElContribution(MatrixDouble &Stiffness) override{
        PanicButton();
    };

    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override;
    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness) override{
        PanicButton();
    };
    void ComputeElContribution(std::vector<VecDouble> &Rhs) override{
        PanicButton();
    };

    void ComputeSpatialDerivatives() override{
        fSubElements[0]->ComputeSpatialDerivatives();
    };

    void ComputeCurrentSpatialDerivatives() override{
        fSubElements[0]->ComputeCurrentSpatialDerivatives();
    };
    void ComputeHighOrderSpatialDerivatives() override{
        fSubElements[0]->ComputeHighOrderSpatialDerivatives();
    };
    double getJacobian() override{
        return fSubElements[0]->getJacobian();
    };


    void interpolateSolution(int &index, VecDouble &u_) override{
        PanicButton();
    };
    void interpolateSolution(VecDouble &phi, VecDouble &u_) override{
        PanicButton();
    };
    void interpolateSolution(){
        for (int i = 0; i < fSubElements.size(); i++){
            fSubElements[i]->interpolateSolution();
        }
    }
    void interpolateSolDerivatives(MatrixDouble &du_dx) override{
        PanicButton();
    }
    void interpolateSolDerivatives(MatrixDouble &dphidx, MatrixDouble &du_dx) override{
        PanicButton();
    }
    void interpolateSolDerivatives() override{
        for (int i = 0; i < fSubElements.size(); i++){
            fSubElements[i]->interpolateSolDerivatives();
        }
    }
    void interpolateSolDTimeDerivatives(VecDouble &du_dt, VecDouble &du_ddt) override{
        PanicButton();
    }
    void interpolateSolDTimeDerivatives() override{
        PanicButton();
    }
    void setIntegPointWeightFunction() override{
        PanicButton();
    }

    int NSides() override{
        return fSubElements[0]->NSides();
    }
    int NSideNodes(int iside) override{
        return fSubElements[0]->NSideNodes(iside);
    }
    int SideNodeLocIndex(int side, int node) override{
        return fSubElements[0]->SideNodeLocIndex(side, node);
    }
    int Dimension() override{
        return fSubElements[0]->Dimension();
    }
    ElementType Type() override{
        return fSubElements[0]->Type();
    }
    int getNumberOfIntegrationPoints() override{
        return fSubElements[0]->getNumberOfIntegrationPoints();
    }
    void ComputeIntPointDistFunction(VecDouble &nodalval) override{
        return fSubElements[0]->ComputeIntPointDistFunction(nodalval);
    }
    double InterpolateVariable(VecDouble &nValues, int point) override{
        return fSubElements[0]->InterpolateVariable(nValues, point);
    }

    void ComputeIntegPointCoordinates() override{
        fSubElements[0]->ComputeIntegPointCoordinates();
    }

};


#endif

