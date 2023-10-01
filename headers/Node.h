//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------NODES-------------------------------------
//------------------------------------------------------------------------------

#ifndef NODE_H
#define NODE_H

// #include <boost/numeric/ublas/io.hpp>

// using namespace boost::numeric;

#include <math.h> 
#include <stdio.h> 
#include <vector>
#include <algorithm>
#include <functional>
#include <array>
#include <iomanip>
#include <cmath>
#include <limits>
#include <ctime>
#include <numeric>
#include <iostream>


#include "DataTypes.h"

/// Defines the node object and stores all nodal variables information
class Node{
private:
    int fNStateVariables;
    VecDouble fSolution;
    VecDouble fPrevSolution;
    VecDouble fCoord;
    VecDouble fPrevCoord;
    VecDouble fInitialCoord;
    VecDouble fDTimeSolution;
    int fIndex, fDimension;
    
    //variables
    VecInt           fConstrainType;        //Constrain direction
    VecDouble        fConstrainValue;       //Nodal prescribed velocity

    VecDouble        fMeshVelocity;           //Nodal mesh velocity

    double           fWeightFunction;         //Nodal Energy Weight Function
    double           fPrevWeightFunction;
    std::vector<int> fInverseIncidence;
    
public:
    Node(VecDouble &coor, int index, int nState){
        fNStateVariables = nState;
        fIndex = index;
        fDimension = coor.size();

        fSolution.resize(fNStateVariables);
        fDTimeSolution.resize(fNStateVariables);
        fPrevSolution.resize(fNStateVariables);
        fCoord.resize(fDimension);
        fPrevCoord.resize(fDimension);
        fInitialCoord.resize(fDimension);
        fPrevSolution.setZero();
        fSolution.setZero();

        fCoord = coor;
        fPrevCoord = coor;
        fInitialCoord = coor;

        fWeightFunction = 0.;   fPrevWeightFunction = 0.;
        
        fInverseIncidence.clear();

        fConstrainValue.resize(fDimension);
        fConstrainType.resize(fDimension);
        fMeshVelocity.resize(fDimension);

        for (int i = 0; i < fDimension; ++i){
            fMeshVelocity[i] = 0.;
            fConstrainType[i] = 0;
            fConstrainValue[i] = 0.;
        }
    }

    int &GetNStateVariables() {return fNStateVariables;}
    void SetSolution(int istate, double sol){
        fPrevSolution[istate] = fSolution[istate];
        fSolution[istate] = sol;
    }
    void IncrementSolution(int istate, double sol){
        fSolution[istate] += sol;
    }
    double GetSolution(int istate){
        return fSolution[istate];
    }

    /// Returns the node coordinate vector
    /// @return node coordinate vector
    VecDouble &getCoordinates() {return fCoord;};

    /// Returns the node coordinate component value
    /// @return node coordinate component value
    double getCoordinateValue(int dir) const {return fCoord[dir];};
    double getPreviousCoordinateValue(int dir) const {return fPrevCoord[dir];};

    /// Returns the node initial coordinate vector
    /// @return node initial coordinate vector
    VecDouble &getInitialCoordinates() {return fInitialCoord;};

    /// Returns the node coordinate vector at the previous time step
    /// @return node coordinate vector at the previous time step
    VecDouble &getPreviousCoordinates() {return fPrevCoord;}

    /// Increment the coordinate vector
    /// @param int direction @param double increment value
    void incrementCoordinate(int dir, double u){fCoord[dir] += u;};

    /// Sets the previous coordinate vector
    /// @param int direction @param double value
    void setPreviousCoordinates(int dir, double u){fPrevCoord[dir] = u;};

    /// Sets the node coordinate vector
    /// @param VecLocD Coordinate
    void setCoordinates(VecDouble &coor){for (int i=0; i<fDimension; i++) fCoord[i] = coor[i];};

    /// Pushs back a term of the inverse incidence, i.e., an element which
    /// contains the node
    /// @param int element
    void pushInverseIncidence(int el) {fInverseIncidence.push_back(el);}

    /// Gets the number of elements which contains the node
    /// @return int number of elements which contains the node
    int getNumberOfElements() const {return fInverseIncidence.size();}

    /// Gets an specific member of the inverse incidence
    /// @param int index @return int element of the inverse incidence
    int getInverseIncidenceElement(int i) const {return fInverseIncidence[i];}
    void clearInverseIncidence(){
        fInverseIncidence.clear();
        fInverseIncidence.shrink_to_fit();
    }
    
    //.........................Mesh Velocity functions..........................
    /// Sets the node mesh velocity
    /// @param double* mesh velocity
    void setMeshVelocity(VecDouble &u){
        for (int i=fDimension; i--; ){
            fMeshVelocity[i] = u[i];          
        };
    };
    void setMeshVelocityComponent(int dir,double u){fMeshVelocity[dir] = u;} ;

    /// Gets the node mesh velocity
    /// @param int direction @return mesh velocity component
    double getMeshVelocity(int dir) const {return fMeshVelocity[dir];}


    //...........................Constrains functions...........................

    /// Sets all node constrains     
    /// @param int direction 
    /// @param int type: 0 - free, 1 - constrained, 2 - glue zone, 
    /// 3 - fluid-structure interface @param double constrain value
    void setConstrains(int dir, int type, double value){
        fConstrainType[dir] = type;
        fConstrainValue[dir] = value;
        fSolution[dir] = value;
    };
    void SetBoundaryCondition(int dir, int type, double value){
        fConstrainType[dir] = type;
        fConstrainValue[dir] = value;
        fSolution[dir] = value;
    };

    void setConstrainValue(int dir, double value){
        fConstrainValue[dir] = value;
    }

    /// Gets node constrain type
    /// @return constrain type
    int getConstrains(int dir) const {return fConstrainType[dir];}

    /// Gets node constrain value
    /// @return constrain value
    double getConstrainValue(int dir) const {return fConstrainValue[dir];}

    //............................Arlequin functions............................
    /// Sets the nodal energy weight function value
    /// @param double weight function value
    void setWeightFunction(double val) {fPrevWeightFunction = fWeightFunction; fWeightFunction = val;};

    /// Gets the nodal energy weight function value
    /// @return weight function value
    double getWeightFunction() const {return fWeightFunction;};
    double getPreviousWeightFunction() const {return fPrevWeightFunction;};
};


#endif

