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
#include<fstream>

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
    VecDouble fDDTimeSolution;
    int64_t fIndex;
    int fDimension;
    VecDouble        fMeshVelocity;           //Nodal mesh velocity
    double           fWeightFunction;         //Nodal Energy Weight Function
    double           fPrevWeightFunction;
    std::vector<int64_t> fInverseIncidence;
    
public:
    Node(VecDouble &coor, int64_t index, int nState = 1);

    void SetNStateVariables(int nstate);

    int &GetNStateVariables() {
        return fNStateVariables;
    }

    void SetSolution(int istate, double sol);

    void IncrementSolution(int istate, double sol);

    VecDouble &Solution() {
        return fSolution;
    }
    
    VecDouble &PrevSolution() {
        return fPrevSolution;
    }

    double GetSolution(int istate){
        return fSolution[istate];
    }

    double GetDSolutionDTime(int istate){
        return fDTimeSolution[istate];
    }
    double GetDSolutionDDTime(int istate){
        return fDDTimeSolution[istate];
    }

    int64_t &Index(){return fIndex;}

    /// Returns the node coordinate vector
    /// @return node coordinate vector
    VecDouble &getCoordinates() {
        return fCoord;
    };

    /// Returns the node coordinate component value
    /// @return node coordinate component value
    double getCoordinateValue(int dir) const {
        return fCoord[dir];
    };
    double getPreviousCoordinateValue(int dir) const {
        return fPrevCoord[dir];
    };

    /// Returns the node initial coordinate vector
    /// @return node initial coordinate vector
    VecDouble &getInitialCoordinates() {
        return fInitialCoord;
    };

    /// Returns the node coordinate vector at the previous time step
    /// @return node coordinate vector at the previous time step
    VecDouble &getPreviousCoordinates() {
        return fPrevCoord;
    }

    /// Increment the coordinate vector
    /// @param int direction @param double increment value
    void incrementCoordinate(int dir, double u);

    /// Sets the previous coordinate vector
    /// @param int direction @param double value
    void setPreviousCoordinates(int dir, double u);

    /// Sets the node coordinate vector
    /// @param VecLocD Coordinate
    void setCoordinates(VecDouble &coor);

    /// Pushs back a term of the inverse incidence, i.e., an element which
    /// contains the node
    /// @param int element
    void pushInverseIncidence(int el);

    /// Gets the number of elements which contains the node
    /// @return int number of elements which contains the node
    int64_t getNumberOfElements() const {
        return fInverseIncidence.size();
    }

    /// Gets an specific member of the inverse incidence
    /// @param int index @return int element of the inverse incidence
    int getInverseIncidenceElement(int i) const {
        return fInverseIncidence[i];
    }

    void clearInverseIncidence();
    
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

