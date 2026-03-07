

//------------------------------------------------------------------------------
//-----------------------------------ELEMENT------------------------------------
//------------------------------------------------------------------------------

#ifndef ELEMENTT_H
#define ELEMENTT_H

#include "Element.h"
#include "Node.h"
#include "PanicButton.h"
#include "IntPointData.h"
// #include "PartitionedQuadrature.hpp"


/// Defines the fluid element object and all the element information
template<class compshape>
class ElementT : public Element{
protected:
    // Integration rule object
    typename compshape::LocIntRule fIntRule;

public:
    ElementT();

    ElementT(int64_t index, GeoElement* gel, CompMesh* mesh, WeakForm *wf);


    double InterpolateVariable(VecDouble &nValues, int point) override;
    void getIntegPointCoordinates();
    
    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    void ComputeSpatialDerivatives() override;
    void ComputeCurrentSpatialDerivatives() override;
    void ComputeHighOrderSpatialDerivatives() override;

    void interpolateSolution(int &index, VecDouble &u_) override;
    void interpolateSolution(VecDouble &phi, VecDouble &u_) override;
    void interpolateSolution() override;
    void interpolateMeshVelocity(int &index, VecDouble &umesh_, VecDouble &umeshPrev_);
    void interpolateSolDerivatives(MatrixDouble &du_dx) override;
    void interpolateSolDerivatives(MatrixDouble &dphidx, MatrixDouble &du_dx) override;
    void interpolateSolDerivatives() override;
    virtual void interpolateSolDTimeDerivatives(VecDouble &du_dt, VecDouble &du_ddt) override;
    virtual void interpolateSolDTimeDerivatives() override;


    void ComputeIntPointDistFunction(VecDouble &nodalval) override;
    /// Gets the element jacobian determinant
    /// @return element jacobinan determinant
    double getJacobian() override {
        // std::cout << "AAA 1 "<< std::endl;
        //Computes the jacobian matrix
        int index = 0;
        fReference->ComputeJacobian(fIntegData);
        // std::cout << "AAA 2 "<< std::endl;        
        //Computes spatial derivatives
        // ComputeSpatialDerivatives(xsi, ainv_, dphi_dx);
        // std::cout << "AAA 3 "<< std::endl;

        return fIntegData.fJacA0;
    };

    // /// Compute and store the drag and lift forces at the element boundary
    // void computeDragAndLiftForces(double &pressureDragForce, double &pressureLiftForce, double &frictionDragForce,
    //                               double &frictionLiftForce, double &dragForce, double &liftForce,
    //                               double &pitchingMoment, double & perimeter);

    // /// Compute and store the boundary forces
    // void getBoundaryLoad(VecDouble &xsi, VecDouble &load);

    

    int Dimension() override {return compshape::Dimension;}  

    

    //.............................Model functions..............................
    //......................Integration Points Information......................
    /// Gets the number of integration points of the special quadrature rule
    /// @retunr number of integration point of the special quadrature rule
    int getNumberOfIntegrationPoints(){return fIntRule.NPoints();};

    void ComputeError(VecDouble &errors) override;
    
    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) override;
    void ComputeElContribution(MatrixDouble &Stiffness) override;
    void ComputeElContribution(VecDouble &Rhs) override;
    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override;
    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness) override;
    void ComputeElContribution(std::vector<VecDouble> &Rhs) override;

    int NSides() override {return compshape::NSides;}
    int NSideNodes(int iside) override {return compshape::NSideNodes(iside);}
    int SideNodeLocIndex(int side, int node) override {return compshape::SideNodeLocIndex(side, node);}

    void setIntegPointWeightFunction() override;

    //Method for creating a copy of the element
    virtual Element *Clone() const;

    

    void Integrate(std::vector<std::string> &varNames, std::map<std::string,VecDouble> &result) override;

    ElementType Type() override{
        return fReference->Type();
    }

};


#endif

