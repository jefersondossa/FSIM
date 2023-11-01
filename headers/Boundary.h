//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------------FLUID BOUNDARY--------------------------------
//------------------------------------------------------------------------------

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "DataTypes.h"
#include "Element.h"

/// Defines the fluid boundary object and its properties
class Boundary : public Element{

public: 
    /// Boundary element constructor
    /// @param BoundConnect boundary element connectivity
    /// @param int boundary element index
    /// @param int boundary element constrain type 
    /// @param int boundary element constrain value @see Node::setConstrains()
    Boundary(VecInt &connec, int index, VecInt &constrain, VecDouble &values, int gr, CompMesh* mesh){
        
        fMesh = mesh;
        fConnect.resize(fMesh->NBdNodes());

        for(int i = 0; i<fMesh->NBdNodes(); i++) fConnect[i] = connec[i];
        fIndex = index;
        group_ = gr;

        constrainType.resize(3);
        constrainValue.resize(3);

        for(int i = 0; i<3; i++){
            constrainType[i] = constrain[i];
            constrainValue[i] = values[i];
        };

        element_ = -1;
        elementSide_ = 1000;
    };     

    /// Returns the boundary element constrain component type
    /// @param int direction @return boundary element constrain component type
    int getConstrain(int dir) override {return constrainType[dir];}

    /// Returns the boundary element constrain component value
    /// @param int direction @return boundary element constrain component value
    double getConstrainValue(int dir) override { return constrainValue[dir];}

    /// Returns the boundary element connectivity
    /// @return boundary element connectivity
    // VecInt &getBoundaryConnectivity(){return connectB_;}
    // void setBoundaryConnectivity(VecInt &connec){for(int i = 0; i<fMesh->NBdNodes(); i++) connectB_[i] = connec[i];}

    /// Sets the boundary element group
    /// @param int boundary element group
    void setBoundaryGroup(int gr) override {group_ = gr;}

    /// Gets the boundary element group
    /// @return boundary element group
    int getBoundaryGroup() override {return group_;};

    /// Sets the fluid element correspondence
    /// @param int fluid element correspondence
    void setElement(int el) override {element_ = el;}

    /// Sets the fluid element correspondence side at the boundary
    /// @param int fluid element correspondence side
    void setElementSide(int el)override {elementSide_ = el;}

    /// Gets the fluid element correspondence
    /// @return fluid element correspondence
    int getElement() override {return element_;}

    /// Gets the fluid element correspondence side at the boundary
    /// @return fluid element correspondence side
    int getElementSide() override {return elementSide_;}

    void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) override{};
    void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) override{};
    void ComputeJacobian(int index) override{};
    void ComputeCurrentJacobian(int index) override{};

    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    void ComputeSpatialDerivatives() override{};
    void ComputeCurrentSpatialDerivatives() override{};
    void interpolateSolution(int &index, VecDouble &u_) override{};
    void interpolateSolution(VecDouble &phi, VecDouble &u_) override{};
    void interpolateSolDerivatives(MatrixDouble &du_dx) override{};

    void setIntersectionParameters(VecDouble &x, VecDouble &X)  override{};
    double getJacobian()  override{ return 0;};
protected:
    
    VecInt       constrainType;     //Element type of constrain
    VecDouble    constrainValue;    //Element constrain value
    int          element_;          //Fluid Element
    int          elementSide_;      //Fluid Element Side
    int          group_;            //Element boundary group
};

#endif