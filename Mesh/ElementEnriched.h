
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

#ifndef ELEMENTENRICHED_H
#define ELEMENTENRICHED_H

#include "ElementT.h"

/// Defines the fluid element object and all the element information
class ElementEnriched : public Element{
private:
    Element* fLocalElement;
    Element* fGlobalElement;
    std::map<int64_t,int64_t> *globalElementCorrespondence;
    std::map<int64_t, MatrixDouble> *globalNodeCorrespondence;
    std::map<int64_t, int64_t> *connectEnrichment;
    
public:
    ElementEnriched() : Element(){
    };

    virtual Element *Clone() const;

    ElementEnriched(int64_t index, Element* localEl, Element* globalEl, CompMesh* mesh, WeakForm *wf);

    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) override;

    void setCorrespondence(std::map<int64_t,int64_t> *elCorresp, std::map<int64_t, MatrixDouble> *nodeCorresp){

        globalElementCorrespondence = elCorresp;
        globalNodeCorrespondence = nodeCorresp;
    };

    void SetEnrichmentData(std::map<int64_t, int64_t> *connectEnrichment_){
        connectEnrichment = connectEnrichment_;
    };

    void ComputeElContribution(VecDouble &Rhs) override{
        PanicButton();
    };
    void ComputeElContribution(MatrixDouble &Stiffness) override{
        PanicButton();
    };

    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override{
        PanicButton();
    };

    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness) override{
        PanicButton();
    };

    void ComputeElContribution(std::vector<VecDouble> &Rhs) override{
        PanicButton();
    };

    void ComputeSpatialDerivatives() override{
        fLocalElement->ComputeSpatialDerivatives();
    };

    void ComputeCurrentSpatialDerivatives() override{
        fLocalElement->ComputeCurrentSpatialDerivatives();
    };

    void ComputeHighOrderSpatialDerivatives() override{
        fLocalElement->ComputeCurrentSpatialDerivatives();
    };

    double getJacobian()override{
        fLocalElement->getJacobian();
    };

    void interpolateSolution(int &index, VecDouble &u_) override{
        PanicButton();
    };

    void interpolateSolution(VecDouble &phi, VecDouble &u_) override{
        PanicButton();
    };

    void interpolateSolution(){
        fLocalElement->interpolateSolution();
    };

    void interpolateSolDerivatives(MatrixDouble &du_dx) override{
        PanicButton();
    };

    void interpolateSolDerivatives(MatrixDouble &dphidx, MatrixDouble &du_dx) override{
        PanicButton();
    };

    void interpolateSolDerivatives() override{
        fLocalElement->interpolateSolDerivatives();
    };

    void interpolateSolDTimeDerivatives(VecDouble &du_dt, VecDouble &du_ddt) override{
        PanicButton();
    };

    void interpolateSolDTimeDerivatives() override{
        PanicButton();
    };

    void setIntegPointWeightFunction() override{
        PanicButton();
    };

    int NSides() override{
        return fLocalElement->NSides();
    };

    int NSideNodes(int iside) override{
        return fLocalElement->NSideNodes(iside);
    };

    int SideNodeLocIndex(int side, int node) override{
        return fLocalElement->SideNodeLocIndex(side, node);
    };

    int Dimension() override{
        return fLocalElement->Dimension();
    };

    ElementType Type() override{
        return fLocalElement->Type();
    };

    int getNumberOfIntegrationPoints() override{
        return fLocalElement->getNumberOfIntegrationPoints();
    };

    double IntegPointWeight(int index) override{
        return fLocalElement->IntegPointWeight(index);
    };

    double IntegPointCoordinate(int index, int dir) override{
        return fLocalElement->IntegPointCoordinate(index, dir);
    };

    void ComputeIntPointDistFunction(VecDouble &nodalval) override{
        return fLocalElement->ComputeIntPointDistFunction(nodalval);
    };

    double InterpolateVariable(VecDouble &nValues, int point) override{
        return fLocalElement->InterpolateVariable(nValues, point);
    };

    void ComputeIntegPointCoordinates() override{
        fLocalElement->ComputeIntegPointCoordinates();
    };

    /*void GetInterfaceData(Element* &localEl, Element* &globalEl) override{
        localEl = fLocalElement;
        globalEl = fGlobalElement;
    };*/
};


#endif

