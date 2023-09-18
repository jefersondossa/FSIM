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
#include "CompMesh.h"

class CompMesh;
/// Defines the fluid boundary object and its properties

class Boundary{

public: 
    /// Boundary element constructor
    /// @param BoundConnect boundary element connectivity
    /// @param int boundary element index
    /// @param int boundary element constrain type 
    /// @param int boundary element constrain value @see Node::setConstrains()
    Boundary(VecInt &connec, int index, VecInt &constrain, VecDouble &values, int gr, CompMesh* mesh){
        
        fMesh = mesh;
        connectB_.resize(fMesh->nBdNodes);

        for(int i = 0; i<fMesh->nBdNodes; i++) connectB_[i] = connec[i];
        index_ = index;
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
    int getConstrain(int dir){return constrainType[dir];}

    /// Returns the boundary element constrain component value
    /// @param int direction @return boundary element constrain component value
    double getConstrainValue(int dir){return constrainValue[dir];}

    /// Returns the boundary element connectivity
    /// @return boundary element connectivity
    VecInt &getBoundaryConnectivity(){return connectB_;}
    void setBoundaryConnectivity(VecInt &connec){for(int i = 0; i<fMesh->nBdNodes; i++) connectB_[i] = connec[i];}

    /// Sets the boundary element group
    /// @param int boundary element group
    void setBoundaryGroup(int gr){group_ = gr;}

    /// Gets the boundary element group
    /// @return boundary element group
    int getBoundaryGroup(){return group_;};

    /// Sets the fluid element correspondence
    /// @param int fluid element correspondence
    void setElement(int el){element_ = el;}

    /// Sets the fluid element correspondence side at the boundary
    /// @param int fluid element correspondence side
    void setElementSide(int el){elementSide_ = el;}

    /// Gets the fluid element correspondence
    /// @return fluid element correspondence
    int getElement(){return element_;}

    /// Gets the fluid element correspondence side at the boundary
    /// @return fluid element correspondence side
    int getElementSide(){return elementSide_;}

private:
    CompMesh* fMesh;
    VecInt       connectB_;         //Boundary element connectivity
    int          index_;            //Boundary element index
    VecInt       constrainType;     //Element type of constrain
    VecDouble    constrainValue;    //Element constrain value
    int          element_;          //Fluid Element
    int          elementSide_;      //Fluid Element Side
    int          group_;            //Element boundary group
};

#endif