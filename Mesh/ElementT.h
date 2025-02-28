
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
#include "IntPointData.h"

/// Defines the fluid element object and all the element information
template<class tshape>
class ElementT : public Element{
protected:
    // Integration rule object
    typename tshape::LocIntRule fIntRule;

public:
    ElementT();

    ElementT(int64_t index, VecInt &connect, CompMesh* mesh, WeakForm *wf);


    double InterpolateVariable(VecDouble &nValues, int point) override;
    void getIntegPointCoordinates();
    //........................Element basic information.........................
    /// Compute and store the spatial jacobian matrix
    /// @param bounded_vector integration point adimensional coordinates
    void ComputeJacobian() override;
    void ComputeCurrentJacobian() override;
    void ComputeJacobianSearch() override;

    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    void ComputeSpatialDerivatives() override;
    void ComputeCurrentSpatialDerivatives() override;
    void getHighOrderSpatialDerivatives(VecDouble &xsi, MatrixDouble &ainv_, MatrixDouble &dphi_dx, MatrixDouble &dDphi_dx);

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
        ComputeJacobian();
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

    

    int Dimension() override {return tshape::Dimension;}  

    //.................Element intersection and correspondence..................
    /// Gets the element intersection parameters 
    /// @param minimum coordinates @param maximum coordinates 
    /// @param minimum inner product @param maximum inner product
    /// @param side lenght
    void setIntersectionParameters(VecDouble &x, VecDouble &X) override ;

    

    //.............................Model functions..............................
    //......................Integration Points Information......................
    /// Gets the number of integration points of the special quadrature rule
    /// @retunr number of integration point of the special quadrature rule
    int getNumberOfIntegrationPoints(){return fIntRule.NPoints();};

    void SetIntPointCoordAndWeight(int index) override;
    VecDouble GetShapeFunction() override;

    int Dimension() const{
        return tshape::Dimension;
    }

    void ComputeError(VecDouble &errors) override;
    
    const int &NElNodes() override {return tshape::NElNodes;}

    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) override;
    void ComputeElContribution(MatrixDouble &Stiffness) override;
    void ComputeElContribution(VecDouble &Rhs) override;
    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override;
    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness) override;
    void ComputeElContribution(std::vector<VecDouble> &Rhs) override;

    int NCornerNodes() override {return tshape::NCornerNodes;}

    void setIntegPointWeightFunction() override;

    //Method for creating a copy of the element
    virtual Element *Clone() const;

    VecDouble NodeCoord(int inode) override{
        VecDouble xnode(3);
        MatrixDouble coords(3,tshape::NElNodes);
        tshape::getCoordinates(coords);
        for (int i = 0; i < 3; i++){
            xnode[i] = coords(i,inode);
        }
        return xnode;
    }

    void Integrate(std::vector<std::string> &varNames, std::map<std::string,VecDouble> &result) override;

    ElementType Type() override{
        return tshape::ElType;
    }

};


#endif

