
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

#ifndef ELEMENTT_H
#define ELEMENTT_H

#include "Element.h"
#include "Node.h"
#include "PanicButton.h"
#include "IntegrationQuadrature.h"
#include "IntegrationQuadrature11.h"
#include "DomainIntegration.h"
#include "IntPointData.h"
// #include "PartitionedQuadrature.hpp"


/// Defines the fluid element object and all the element information
template<class tshape>
class ElementT : public Element{
protected:
    // Integration rule object
    typename tshape::LocIntRule fIntRule;

public:
    ElementT() : Element(){
        fIntRule.SetOrder(tshape::Order+tshape::Order);
    };

    /// fluid element constructor
    /// @param int element index @param Connectivity element connectivity
    /// @param vector<Nodes> 
    ElementT(int64_t index, VecInt &connect, CompMesh* mesh) : Element(){
        
        fMesh = mesh;
        fConnect.resize(tshape::NElNodes);
        fIndex = index;
        for (int i = tshape::NElNodes; i--; ) fConnect[i] = connect[i];
        DEG = fMesh->GetDefaultOrder();

        FSIInterface = false;
        fNeighborElements.clear();

        fIntRule.SetOrder(tshape::Order+tshape::Order);
        fIntegData.fWeightFunction.resize(fIntRule.NPoints());
        fIntegData.fDistFunction.resize(fIntRule.NPoints());
        fIntegData.fPrevWeightFunction.resize(fIntRule.NPoints());

        fIntegData.fWeightFunction.fill(1.);
        fIntegData.fPrevWeightFunction.fill(1.);
    
        getIntegPointCoordinates();

    };

    ElementT(int64_t index, VecInt &connect, CompMesh* mesh, WeakForm *wf) : Element(){
        
        fMesh = mesh;
        fConnect.resize(tshape::NElNodes);
        fIndex = index;
        for (int i = tshape::NElNodes; i--; ) fConnect[i] = connect[i];
        DEG = fMesh->GetDefaultOrder();
        fWeakForm = wf;
        if (fWeakForm) nLocDOF = tshape::NElNodes * fWeakForm->NState(); 

        FSIInterface = false;
        fNeighborElements.clear();
        int increase = 0;
        if(wf->GetExactSolution()) increase = 2;

        fIntRule.SetOrder(tshape::Order+tshape::Order+increase);

        fIntegData.fWeightFunction.resize(fIntRule.NPoints());
        fIntegData.fDistFunction.resize(fIntRule.NPoints());
        fIntegData.fPrevWeightFunction.resize(fIntRule.NPoints());

        fIntegData.fWeightFunction.fill(1.);
        fIntegData.fPrevWeightFunction.fill(1.);
    
        getIntegPointCoordinates();

    };


    double InterpolateVariable(VecDouble &nValues, int point) override;
    void getIntegPointCoordinates();
    //........................Element basic information.........................
    /// Clear all element variables
    void clearVariables();

    /// Compute and store the spatial jacobian matrix
    /// @param bounded_vector integration point adimensional coordinates
    void ComputeJacobian() override;
    void ComputeCurrentJacobian() override;

    
    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    void ComputeSpatialDerivatives() override;
    void ComputeCurrentSpatialDerivatives() override;
    void getHighOrderSpatialDerivatives(VecDouble &xsi, MatrixDouble &ainv_, MatrixDouble &dphi_dx, MatrixDouble &dDphi_dx);

    void interpolateSolution(int &index, VecDouble &u_) override;
    void interpolateSolution(VecDouble &phi, VecDouble &u_) override;
    void interpolateSolution();
    void interpolateMeshVelocity(int &index, VecDouble &umesh_, VecDouble &umeshPrev_);
    void interpolateSolDerivatives(MatrixDouble &du_dx) override;
    void interpolateSolDerivatives(MatrixDouble &dphidx, MatrixDouble &du_dx) override;
    void interpolateSolDerivatives();


    void ComputeIntPointDistFunction(VecDouble &nodalval) override;
    /// Gets the element jacobian determinant
    /// @return element jacobinan determinant
    double getJacobian() override {
        // std::cout << "AAA 1 "<< std::endl;
        //Computes the jacobian matrix
        int index = 0;
        ComputeJacobian();
        // std::cout << "AAA 2 "<< std::endl;        
        //Computes spatial derivatives
        // ComputeSpatialDerivatives(xsi, ainv_, dphi_dx);
        // std::cout << "AAA 3 "<< std::endl;

        return fIntegData.fJacA0;
    };

    /// Compute and store the drag and lift forces at the element boundary
    void computeDragAndLiftForces(double &pressureDragForce, double &pressureLiftForce, double &frictionDragForce,
                                  double &frictionLiftForce, double &dragForce, double &liftForce,
                                  double &pitchingMoment, double & perimeter);

    /// Compute and store the boundary forces
    void getBoundaryLoad(VecDouble &xsi, VecDouble &load);

    

    int Dimension() override {return tshape::Dimension;}  

    

    /// Gets the boundary connectivity for boundary integration
    /// @param int* boundary connectivity
    void getBoundaryNodes(int *nodesb_);

    //.................Element intersection and correspondence..................
    /// Gets the element intersection parameters 
    /// @param minimum coordinates @param maximum coordinates 
    /// @param minimum inner product @param maximum inner product
    /// @param side lenght
    void setIntersectionParameters(VecDouble &x, VecDouble &X) override ;

    

    //.............................Model functions..............................
    /// Sets if the element belongs to the fluid structure interface
    
    bool &getFSIInterface(){return FSIInterface;};

    //......................Integration Points Information......................
    /// Gets the number of integration points of the special quadrature rule
    /// @retunr number of integration point of the special quadrature rule
    int getNumberOfIntegrationPoints(){return fIntRule.NPoints();};

    int Dimension() const{
        return tshape::Dimension;
    }

    
    const int &NElNodes() override {return tshape::NElNodes;}

    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) override;
    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override;

    int NCornerNodes() override {return tshape::NCornerNodes;}

    void setIntegPointWeightFunction() override;

    int getElement()override{return 0;};
    int getElementSide()override {return 0;};
    void setElement(int el)override{ };
    void setElementSide(int el)override {};

    //Method for creating a copy of the element
    virtual Element *Clone() const;
};


#endif

