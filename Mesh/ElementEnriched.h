
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
#include "MixedCompMesh.h"

/// Defines the fluid element object and all the element information
class ElementEnriched : public Element{
private:
    Element* fLocalElement;
    Element* fGlobalElement;
    std::map<int,int> *globalElementCorrespondence;
    std::map<int, MatrixDouble> *globalNodeCorrespondence;
    std::map<int, int> *connectEnrichment;
    
public:
    ElementEnriched() : Element(){
    };

    virtual Element *Clone() const;

    ElementEnriched(int index, Element* localEl, Element* globalEl, CompMesh* mesh, WeakForm *wf);

    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) override;

    void setCorrespondence(std::map<int,int> *elCorresp, std::map<int, MatrixDouble> *nodeCorresp){

        globalElementCorrespondence = elCorresp;
        globalNodeCorrespondence = nodeCorresp;
    };

    void SetEnrichmentData(std::map<int, int> *connectEnrichment_){
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

    REAL getJacobian()override{
        fLocalElement->getJacobian();
    };

    void interpolateSolution(int &index, VecDouble &u_) override{
        PanicButton();
    };

    void interpolateSolution(VecDouble &phi, VecDouble &u_) override{
        PanicButton();
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

    REAL IntegPointWeight(int index) override{
        return fLocalElement->IntegPointWeight(index);
    };

    REAL IntegPointCoordinate(int index, int dir) override{
        return fLocalElement->IntegPointCoordinate(index, dir);
    };

    void ComputeIntPointDistFunction(VecDouble &nodalval) override{
        return fLocalElement->ComputeIntPointDistFunction(nodalval);
    };

    REAL InterpolateVariable(VecDouble &nValues, int point) override{
        return fLocalElement->InterpolateVariable(nValues, point);
    };

    void ComputeIntegPointCoordinates() override{
        fLocalElement->ComputeIntegPointCoordinates();
    };

    IntPointData &IntegrationData()override{
        return fLocalElement->IntegrationData();
    }

    void interpolateSolution(){
        // fLocalElement->IntegrationData().fAdimCoord = fIntegData.fAdimCoord;
        // fLocalElement->IntegrationData().fWeight = fIntegData.fWeight;
        // fLocalElement->IntegrationData().fA0 = this->fIntegData.fA0;
        // fLocalElement->IntegrationData().fA0Inv = this->fIntegData.fA0Inv;
        // fLocalElement->IntegrationData().fAxes0 = this->fIntegData.fAxes0;
        // fLocalElement->IntegrationData().fJacA0 = this->fIntegData.fJacA0;
        // fLocalElement->IntegrationData().fX = this->fIntegData.fX;

        // //Computes spatial derivatives
        // fLocalElement->ComputeSpatialDerivatives();
        // fLocalElement->interpolateSolution();
        // fLocalElement->interpolateSolDerivatives();

        std::cout << "Para interpolar a solução global enriquecida é necessário encontrar também as correspondências dos nós da malha local\nInserir o valor correto aqui:";
        
        for (int k = 0; k < fLocalElement->Dimension(); k++){   
            fGlobalElement->IntegrationData().fAdimCoord[k] = globalNodeCorrespondence->at(fLocalElement->Index())(fLocalElement->IntegrationData().fIndex,k);
        }
        // fGlobalElement->IntegrationData().fAdimCoord = this->fIntegData.fAdimCoord;
        fGlobalElement->IntegrationData().fWeight = this->fIntegData.fWeight;
        fGlobalElement->Reference()->ComputeJacobian(fGlobalElement->IntegrationData());

        //Computes spatial derivatives
        fGlobalElement->ComputeSpatialDerivatives();
        fGlobalElement->interpolateSolution();
        fGlobalElement->interpolateSolDerivatives();
        
        
        // PanicButton();
    }

    void Solution(int var, VecDouble &Sol) override{
        
        fWeakForm->Solution(fGlobalElement->IntegrationData(),var,Sol);
    };

    int NShapeFunctions() override{
        return fGlobalElement->NShapeFunctions();
    };

   void setConnectivity(int mesh, VecInt connect){
        fConnect.resize(connect.size());
        MixedCompMesh* mixedMesh = dynamic_cast<MixedCompMesh*>(fMesh);
        for (size_t i = 0; i < connect.size(); i++){
            fConnect[i] = mixedMesh->MeshVector()[mesh]->ConnectVec()[connect[i]];
        }
    };

    void setConnectivity(VecInt connect) override{
        for (int iconnect = 0; iconnect < connect.size(); iconnect++){
            fConnect[iconnect] = fMesh->ConnectVec()[connect[iconnect]];
        } 
    }

    void setConnectivity(std::vector<Connect *> connect)override{fConnect = connect;};


    void AccountForEnrichment(MatrixDouble &Stiffness, VecDouble &Rhs);
    void AccountForEnrichmentMixed(MatrixDouble &Stiffness, VecDouble &Rhs);

};


#endif

