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

    std::vector<Element *> &SubElements(){return fSubElements;};

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
        this->fIntegData.fSol = fSubElements[0]->IntegrationData().fSol;
    }

    void Solution(int var, VecDouble &Sol) override{
        for (int i = 0; i < fSubElements.size(); i++){
            fSubElements[i]->IntegrationData().fAdimCoord = fIntegData.fAdimCoord;
            fSubElements[i]->IntegrationData().fWeight = fIntegData.fWeight;
            fSubElements[i]->IntegrationData().fA0 = this->fIntegData.fA0;
            fSubElements[i]->IntegrationData().fA0Inv = this->fIntegData.fA0Inv;
            fSubElements[i]->IntegrationData().fAxes0 = this->fIntegData.fAxes0;
            fSubElements[i]->IntegrationData().fJacA0 = this->fIntegData.fJacA0;
            fSubElements[i]->IntegrationData().fX = this->fIntegData.fX;

            //Computes spatial derivatives
            fSubElements[i]->ComputeSpatialDerivatives();
            fSubElements[i]->interpolateSolution();
            fSubElements[i]->interpolateSolDerivatives();
        }
        std::vector<IntPointData *> data(fSubElements.size());
        for (int i = 0; i < fSubElements.size(); i++){
            data[i] = &fSubElements[i]->IntegrationData();
        }
        fWeakForm->Solution(data,var,Sol);
    };

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
        this->fIntegData.fDSolDx = fSubElements[0]->IntegrationData().fDSolDx;
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
    double IntegPointWeight(int index) override{
        return fSubElements[0]->IntegPointWeight(index);
    }
    double IntegPointCoordinate(int index, int dir) override{
        return fSubElements[0]->IntegPointCoordinate(index, dir);
    }
    void ComputeIntPointDistFunction(VecDouble &nodalval) override{
        return fSubElements[0]->ComputeIntPointDistFunction(nodalval);
    }
    double InterpolateVariable(VecDouble &nValues, int point) override{
        return fSubElements[0]->InterpolateVariable(nValues, point);
    }

    void ComputeIntegPointCoordinates() override{
        fSubElements[0]->ComputeIntegPointCoordinates();
        this->fIntPointCoordinates = fSubElements[0]->getIntPointCoordinates();
    }

    int NShapeFunctions() override{
        int nshape = 0;
        for (int i = 0; i < fSubElements.size(); i++){
            nshape += fSubElements[i]->NShapeFunctions();
        }
        return nshape;
    };

};


#endif

